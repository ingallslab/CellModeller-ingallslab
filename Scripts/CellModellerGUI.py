# CellModeller4
#
#
# GUI
#
# Tim Rudge
# Jan 2011

import sys
from importlib.resources import files

from PyQt5 import uic
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication


def setup_ui():
    # The Qt application
    qapp = QApplication([])

    # The UI
    fileref = files("CellModeller.GUI").joinpath("PyGLGUI.ui")
    with fileref.open("rb") as uifile:
        ui = uic.loadUi(uifile)

        if ui is not None:
            ui.show()
            ui.raise_()

            label = ui.label
            label.setTextFormat(Qt.TextFormat.RichText)
            label.setAlignment(Qt.AlignmentFlag.AlignJustify)

            viewer_widget = ui.PyGLCMViewer
            # Load a model if specified
            if len(sys.argv) > 1:
                viewer_widget.loadModelFile(sys.argv[1])

            # Launch app main loop
            sys.exit(qapp.exec_())


if __name__ == "__main__":
    setup_ui()
