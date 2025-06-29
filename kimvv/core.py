import os
import pathlib

import kim_edn


class KIMVVTestDriver:
    @property
    def kimspec(self):
        mypath = pathlib.Path(__file__).parent.resolve()
        myname = self.__class__.__name__
        return kim_edn.load(os.path.join(mypath, myname, "kimspec.edn"))
