from .Mygen import GEMOptions

class ConfOpt:
    def __init__(self, **kwargs):
        self._opt = GEMOptions()  # private internal config object
        self.set(**kwargs)

    def set(self, **kwargs):
        for k, v in kwargs.items():
            if v is None:
                continue
            if hasattr(self._opt, k):
                setattr(self._opt, k, v)
            else:
                raise ValueError(f"Invalid GEMOptions field: '{k}'")

        # if self._opt.threads > 0 and self._opt.num_chunks > 0:
        #     if self._opt.threads > self._opt.num_chunks:
        #         raise ValueError(
        #             f"Threads {self._opt.threads} cannot be greater than num_chunks {self._opt.num_chunks}"
        #         )
    # def __getattr__(self, name):
    #     # only called if attribute not found in ConfOpt
    #     return getattr(self._opt, name)


    def get(self):
        return self._opt
    
    @property
    def outfile(self):
        return self._opt.outfile
        
    @property
    def chunk_size(self):
        return self._opt.num_chunks
    
    @property
    def threads(self):
        return self._opt.threads
        
    @property
    def stream_snps(self):
        return self._opt.stream_snps
