"""
Maybe this feature can be reimplemented as a CMModule?
Unsure will have to think what needs to be saved to a file...
"""

import inspect
import os
import sys
import time

from CellModeller.Modules.CMModule import CMModule


class PicklesWriter(CMModule):
    def __init__(self, moduleName: str, moduleStr: str = "", dirName: str = "") -> None:
        if dirName:
            outputFileRoot = dirName
        else:
            outputFileRoot = moduleName + "-"
            outputFileRoot += time.strftime("%y-%m-%d-%H-%M", time.localtime())
        self.outputDirPath = os.path.join("data", outputFileRoot)

        # Add a number to end of dir name if it already exists
        label = 2
        while os.path.exists(self.outputDirPath):
            if label > 2:
                self.outputDirPath = self.outputDirPath[:-2] + "_" + str(label)
            else:
                self.outputDirPath = self.outputDirPath + "_" + str(label)
            label += 1
        os.makedirs(self.outputDirPath)

        # write a copy of the model into the dir (for reference),
        # this goes in the pickle too (and gets loaded when a pickle is loaded)
        if not moduleStr:
            moduleStr = inspect.getsource(sys.modules[moduleName])
        open(os.path.join(self.outputDirPath, moduleName), "w").write(moduleStr)

    def on_sim_ready(self) -> None:
        return None

    def sim_step(self, dt, cell_states) -> None:
        """
        This is where the saving is implemented?
        """
        del dt, cell_states
