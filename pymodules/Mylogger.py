import sys, os, threading, atexit

class FDTee:
    """
    Redirect process stdout+stderr to a log file.
    Optionally also mirror to terminal.
    Captures Python prints AND C++ std::cout/std::cerr.
    """
    def __init__(self, log_path: str, truncate: bool = False, tee_to_terminal: bool = False):
        flags = os.O_WRONLY | os.O_CREAT | (os.O_TRUNC if truncate else os.O_APPEND)
        self.log_fd = os.open(log_path, flags, 0o644)

        self.tee_to_terminal = tee_to_terminal

        # Save originals so we can restore later
        self.orig_out = os.dup(1)
        self.orig_err = os.dup(2)

        # Create pipe and redirect stdout/stderr to it
        self.r_fd, self.w_fd = os.pipe()
        os.dup2(self.w_fd, 1)
        os.dup2(self.w_fd, 2)
        os.close(self.w_fd)

        self._closed = False
        self._lock = threading.Lock()

        self.t = threading.Thread(target=self._pump, daemon=True)
        self.t.start()

        # Make Python wrappers flush aggressively
        try:
            sys.stdout.reconfigure(line_buffering=True, write_through=True)
            sys.stderr.reconfigure(line_buffering=True, write_through=True)
        except Exception:
            pass

        atexit.register(self.close)

    def _safe_write(self, fd: int, data: bytes) -> bool:
        try:
            os.write(fd, data)
            return True
        except OSError:
            return False

    def _pump(self):
        while True:
            try:
                data = os.read(self.r_fd, 65536)  
            except OSError:
                break
            if not data:
                break

            # Always write to log
            self._safe_write(self.log_fd, data)

            # Only write to terminal if requested
            if self.tee_to_terminal:
                self._safe_write(self.orig_out, data)

        try:
            os.close(self.r_fd)
        except OSError:
            pass

    def close(self):
        with self._lock:
            if self._closed:
                return
            self._closed = True

            try:
                sys.stdout.flush()
                sys.stderr.flush()
            except Exception:
                pass

            # Restore stdout/stderr (closes pipe write end)
            try:
                os.dup2(self.orig_out, 1)
                os.dup2(self.orig_err, 2)
            except OSError:
                pass

        if self.t.is_alive():
            self.t.join(timeout=2.0)

        for fd in (self.orig_out, self.orig_err, self.log_fd):
            try:
                os.close(fd)
            except OSError:
                pass
