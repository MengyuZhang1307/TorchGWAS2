import os, threading, atexit

class FDTee:
    """
    Tee process-level stdout+stderr to BOTH terminal and log file.
    Captures Python prints AND C++ std::cout/std::cerr.
    """
    def __init__(self, log_path: str, truncate: bool = False):
        flags = os.O_WRONLY | os.O_CREAT | (os.O_TRUNC if truncate else os.O_APPEND)
        self.log_fd = os.open(log_path, flags, 0o644)

        self.orig_out = os.dup(1)  # terminal/pipeline stdout
        self.orig_err = os.dup(2)  # terminal/pipeline stderr

        self.r_fd, self.w_fd = os.pipe()

        # Redirect process stdout/stderr -> pipe
        os.dup2(self.w_fd, 1)
        os.dup2(self.w_fd, 2)
        os.close(self.w_fd)

        self._closed = False
        self._lock = threading.Lock()

        self.t = threading.Thread(target=self._pump, daemon=True)
        self.t.start()

        atexit.register(self.close)

    def _safe_write(self, fd: int, data: bytes) -> bool:
        try:
            os.write(fd, data)
            return True
        except OSError:
            return False

    def _pump(self):
        write_term = True
        write_log = True
        while True:
            try:
                data = os.read(self.r_fd, 4096)
            except OSError:
                break

            if not data:
                break

            # terminal
            if write_term:
                ok = self._safe_write(self.orig_out, data)
                if not ok:
                    write_term = False  # stop trying terminal

            # log file
            if write_log:
                ok = self._safe_write(self.log_fd, data)
                if not ok:
                    write_log = False  # stop trying log

            # if both destinations failed, still keep draining (discard output)
            # so the main process never blocks.

        try:
            os.close(self.r_fd)
        except OSError:
            pass

    def close(self):
        with self._lock:
            if self._closed:
                return
            self._closed = True

            # restore stdout/stderr
            try:
                os.dup2(self.orig_out, 1)
                os.dup2(self.orig_err, 2)
            except OSError:
                pass

            for fd in (self.orig_out, self.orig_err, self.log_fd):
                try:
                    os.close(fd)
                except OSError:
                    pass
