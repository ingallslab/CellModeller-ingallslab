"""
@brief Namespace for main CellModeller GUI class.
@details Contains constants relevant to GUI window and rendering.
"""

import sys
import time
from collections.abc import Callable
from typing import Any, Optional

import OpenGL.GL as gl
import pyopencl as cl
from PyQt5.QtCore import (
    QObject,
    QPoint,
    Qt,
    QThread,
    QTimer,
    pyqtBoundSignal,
    pyqtSignal,
    pyqtSlot,
)
from PyQt5.QtGui import (
    QKeyEvent,
    QMatrix4x4,
    QMouseEvent,
    QSurfaceFormat,
    QVector3D,
    QWheelEvent,
)
from PyQt5.QtWidgets import (
    QFileDialog,
    QInputDialog,
    QMessageBox,
    QOpenGLWidget,
    QWidget,
)

from CellModeller.GUI.Renderers.BaseRenderers.GridRenderer import GridRenderer
from CellModeller.GUI.Renderers.Renderer import Renderer, WriteFunc
from CellModeller.Simulator import ModuleTypeError, Simulator

## @brief Time step for simulation in seconds.
DT = 0.05
## @brief Attempt to perform render loop this often per second.
MAX_FRAME_RATE = 60

## @name UI constants
# Change how users interact with GUI window.
# @{

## @brief Scale mouse coordinates to world units.
MOUSE_SCALE_FACTOR = 1 / 200
## @brief Scale mouse wheel to world units.
MOUSE_SCROLL_FACTOR = 1 / 100
## @brief (Translational) in world units/second.
KB_MVMT_SPD = 10
## @brief In radians/second ...uhhh probably.
KB_ROT_SPD = 2
## @brief Speed modifier when holding key (default: Shift).
KB_SLOW_SPD = 0.25
## @brief Relative modifier (applies to both keyboard and mouse).
HORZ_ROT_SPD = 0.5
## @brief see above.
VERT_ROT_SPD = 0.5
## @brief see above.
TILT_ROT_SPD = 0.5
## @brief see above.
HORZ_MVMT_SPD = 5
## @brief see above.
VERT_MVMT_SPD = 5
## @brief see above (depth axis).
Z_MVMT_SPD = 5
## @brief For projection matrix (does this actually work?).
FOV = 45
## @brief Closest scroll zoom distance.
NEAR_DIST = 1
## @brief Farthest scroll zoom distance.
FAR_DIST = 1000

## @}


class PyGLCMViewer(QOpenGLWidget):
    """
    @brief GL widget for CellModeller GUI.
    """

    def keybindings(self) -> None:
        """
        @brief Default bindings for keyboard control.
        @details
        Key     | Function
        --------|----------
        R       | Reset view to initial position.
        G       | Toggle grid visibility.
        P       | Snap view rotation to closest axis.
        W/A/S/D | Move up/left/down/right relative to view.
        Arrows  | ^
        E/Q     | Move closer/further relative to view.
        PgUp/PgDown| ^
        I/J/K/L | Rotate view up/left/down/right.
        U/O     | Tilt view left/right.
        """
        self.single_key_bindings: dict[int, Callable[[], Any]] = {
            Qt.Key.Key_R: lambda: self._reset_view(),
            Qt.Key.Key_G: lambda: (
                self.renderers.update({"Grid": (self.grid, 1)})
                if "Grid" not in self.renderers
                else self.renderers.pop("Grid")
            ),
            Qt.Key.Key_P: lambda: self._snap_up(),
        }
        self.press_key_bindings: dict[int, Callable[[float], None]] = {
            Qt.Key.Key_W: lambda dt: self._move_view(0, 0, dt * KB_MVMT_SPD),
            Qt.Key.Key_Up: lambda dt: self._move_view(0, 0, dt * KB_MVMT_SPD),
            Qt.Key.Key_A: lambda dt: self._move_view(dt * KB_MVMT_SPD, 0, 0),
            Qt.Key.Key_Left: lambda dt: self._move_view(dt * KB_MVMT_SPD, 0, 0),
            Qt.Key.Key_S: lambda dt: self._move_view(0, 0, -dt * KB_MVMT_SPD),
            Qt.Key.Key_Down: lambda dt: self._move_view(0, 0, -dt * KB_MVMT_SPD),
            Qt.Key.Key_D: lambda dt: self._move_view(-dt * KB_MVMT_SPD, 0, 0),
            Qt.Key.Key_Right: lambda dt: self._move_view(-dt * KB_MVMT_SPD, 0, 0),
            Qt.Key.Key_E: lambda dt: self._move_view(0, dt * KB_MVMT_SPD, 0),
            Qt.Key.Key_PageUp: lambda dt: self._move_view(0, dt * KB_MVMT_SPD, 0),
            Qt.Key.Key_Q: lambda dt: self._move_view(0, -dt * KB_MVMT_SPD, 0),
            Qt.Key.Key_PageDown: lambda dt: self._move_view(0, -dt * KB_MVMT_SPD, 0),
            Qt.Key.Key_I: lambda dt: self._rotate_view(0, -dt * KB_ROT_SPD, 0),
            Qt.Key.Key_J: lambda dt: self._rotate_view(-dt * KB_ROT_SPD, 0, 0),
            Qt.Key.Key_K: lambda dt: self._rotate_view(0, dt * KB_ROT_SPD, 0),
            Qt.Key.Key_L: lambda dt: self._rotate_view(dt * KB_ROT_SPD, 0, 0),
            Qt.Key.Key_U: lambda dt: self._rotate_view(0, 0, dt * KB_ROT_SPD),
            Qt.Key.Key_O: lambda dt: self._rotate_view(0, 0, -dt * KB_ROT_SPD),
            Qt.Key.Key_Shift: lambda _: None,
        }

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)

        # Setup OpenGL context profile
        format = QSurfaceFormat()
        format.setDepthBufferSize(24)
        format.setStencilBufferSize(8)
        self.setFormat(format)

        # accept keyboard focus
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

        # Starting values
        self.sim: Simulator | None = None
        self.app_state = {"run": False}
        self.clPlatformNum = 0
        self.clDeviceNum = 0

        # smooth keyboard movement
        self.keys_pressed: dict[int, Callable[[float], None]] = dict()
        self.last_time = time.time()
        # incremental update time (consumed by sim.step)
        self.delta = 0.0
        # average times for fps counter
        self.lastnframes: list[float] = []
        # For click and drag behaviour
        self.last_mouse_pt = QPoint()

        # Initialize animation loop
        self.event_loop = QTimer()
        self.event_loop.setTimerType(Qt.TimerType.PreciseTimer)
        self.event_loop.timeout.connect(self._animate)
        self.event_loop.start(int(1000 / MAX_FRAME_RATE))

        # simulation update thread
        self.sim_thread: QThread | None = None

        # include renderer priority
        # grid needs to be drawn last to avoid early discarding in depth test!
        self.grid = GridRenderer()
        self.default_renderers: dict[str, tuple[Renderer, int]] = {
            "Grid": (self.grid, 1)
        }
        self.renderers: dict[str, tuple[Renderer, int]] = self.default_renderers
        self.buffers_to_write: dict[str, list[WriteFunc]] = {}

        # init coordinate transformation matrices
        self._reset_view()

    ## @name User interaction
    # Handle user inputs - keyboard and mouse.
    # @{

    def keyPressEvent(self, a0: Optional[QKeyEvent]) -> None:
        """
        @brief from @c QtWidgets.QOpenGLWidget.
        """
        if a0 is not None:
            if a0.key() in self.single_key_bindings:
                self.single_key_bindings[a0.key()]()
            elif a0.key() in self.press_key_bindings:
                self.keys_pressed[a0.key()] = self.press_key_bindings[a0.key()]
            a0.accept()

    def keyReleaseEvent(self, a0: Optional[QKeyEvent]) -> None:
        """
        @brief from @c QtWidgets.QOpenGLWidget.
        """
        if a0 is not None:
            if a0.key() in self.keys_pressed:
                self.keys_pressed.pop(a0.key())

    def _handle_keys(self, dt: float) -> None:
        if Qt.Key.Key_Shift in self.keys_pressed:
            dt *= KB_SLOW_SPD
        for handler in self.keys_pressed.values():
            handler(dt)

    def _get_fps(self, dt: float) -> None:
        n = 10
        self.lastnframes.insert(0, dt)
        if len(self.lastnframes) >= n:
            fps = sum(self.lastnframes)
            fps /= n
            if fps != 0:
                fps = int(1 / fps)
                self.selectedCell.emit(f"FPS: {fps}")
            self.lastnframes.pop()

    def wheelEvent(self, a0: Optional[QWheelEvent]) -> None:
        """
        @brief from @c QtWidgets.QOpenGLWidget.
        """
        if a0 is not None:
            degrees = a0.angleDelta()
            if not degrees.isNull():
                self._zoom(degrees.y() * MOUSE_SCROLL_FACTOR)
            a0.accept()

    def mousePressEvent(self, a0: Optional[QMouseEvent]) -> None:
        """
        @brief from @c QtWidgets.QOpenGLWidget.
        """
        if a0 is not None:
            self.last_mouse_pt = a0.pos()
            a0.accept()

    def mouseMoveEvent(self, a0: Optional[QMouseEvent]) -> None:
        """
        @brief from @c QtWidgets.QOpenGLWidget.
        """
        if a0 is not None:
            delta = a0.pos() - self.last_mouse_pt
            if a0.buttons() & Qt.MouseButton.LeftButton == Qt.MouseButton.LeftButton:
                self._rotate_view(
                    delta.x() * MOUSE_SCALE_FACTOR, delta.y() * MOUSE_SCALE_FACTOR, 0
                )
            elif (
                a0.buttons() & Qt.MouseButton.RightButton == Qt.MouseButton.RightButton
            ):
                self._move_view(
                    delta.x() * MOUSE_SCALE_FACTOR, 0, delta.y() * MOUSE_SCALE_FACTOR
                )
            self.last_mouse_pt = a0.pos()
            a0.accept()

    ## @}

    ## @name Camera functions
    # Set matrices for OpenGL.
    # @{

    def _reset_view(self) -> None:
        """
        @brief Set to default.
        """
        self.view_pos = QVector3D(0, 0, 32)
        self.view_dir = QVector3D(0, 0, 0)
        self.view_up = QVector3D(0, 1, 0)
        dist = self.view_pos - self.view_dir
        up = self.view_up.normalized()
        right = QVector3D.normal(dist, up)
        self.view_up = QVector3D.normal(right, dist)

        self.view = QMatrix4x4()
        self.view.lookAt(self.view_pos, self.view_dir, self.view_up)
        self.view_inv = self.view.inverted()[0]
        self.projection = QMatrix4x4()
        self.projection.perspective(
            FOV, self.width() / self.height(), NEAR_DIST, FAR_DIST
        )
        self.projection_inv = self.projection.inverted()[0]

    def _snap_up(self) -> None:
        """
        @brief Project 'up' vector to closest axis.
        """
        x = self.view_up.x()
        y = self.view_up.y()
        z = self.view_up.z()
        coords = [x, y, z]
        abs_coords = list(map(abs, coords))
        i = abs_coords.index(max(abs_coords))
        new_up = [0, 0, 0]
        new_up[i] = 1 if coords[i] >= 0 else -1
        self.view_up = QVector3D(*new_up)

        self.view.setToIdentity()
        self.view.lookAt(self.view_pos, self.view_dir, self.view_up)
        self.view_inv = self.view.inverted()[0]

    def _zoom(self, d: float) -> None:
        """
        @brief Change camera distance from center of view.
        """
        dist = self.view_pos - self.view_dir
        old_length = dist.length()
        new_length = old_length + d
        if new_length <= NEAR_DIST or new_length >= FAR_DIST:
            return

        dist *= new_length / old_length
        self.view_pos = self.view_dir + dist

        self.view.setToIdentity()
        self.view.lookAt(self.view_pos, self.view_dir, self.view_up)
        self.view_inv = self.view.inverted()[0]

    def _rotate_view(self, dx: float, dy: float, dz: float) -> None:
        """
        @brief Change orientation by specified amount in each coordinate.
        """
        dist = self.view_pos - self.view_dir
        old_length = dist.length()
        up = self.view_up.normalized()
        right = QVector3D.normal(dist, up)

        self.view_up = QVector3D.normal(right, dist)
        self.view_up += right * dz * TILT_ROT_SPD
        self.view_up.normalize()

        right *= dx * HORZ_ROT_SPD * old_length
        up *= dy * VERT_ROT_SPD * old_length
        dist += up + right
        self.view_dir = self.view_pos - dist.normalized() * old_length

        self.view.setToIdentity()
        self.view.lookAt(self.view_pos, self.view_dir, self.view_up)
        self.view_inv = self.view.inverted()[0]

    def _move_view(self, dx: float, dy: float, dz: float) -> None:
        """
        @brief Translate camera by specified amount in each coordinate.
        """
        dist = self.view_pos - self.view_dir
        up = self.view_up.normalized()
        right = QVector3D.normal(dist, up)
        forward = QVector3D.normal(right, up)

        up *= dz * Z_MVMT_SPD
        right *= dx * HORZ_MVMT_SPD
        forward *= dy * VERT_MVMT_SPD
        self.view_pos += up + right + forward
        self.view_dir += up + right + forward

        self.view.setToIdentity()
        self.view.lookAt(self.view_pos, self.view_dir, self.view_up)
        self.view_inv = self.view.inverted()[0]

    ## @}

    ## @name OpenGL functions
    # QOpenGLWidget virtual functions to manage OpenGL context
    # @{

    def initializeGL(self) -> None:
        """
        @brief from @c QtWidgets.QOpenGLWidget.
        """
        self.initialize_renderers()

    def initialize_renderers(self) -> None:
        # automatically called in initializeGL() but necessary otherwise
        self.makeCurrent()

        if self.sim is not None:
            max_cells = self.sim.user_module.max_cells
        else:
            max_cells = 1

        self.renderers = dict(sorted(self.renderers.items(), key=lambda x: x[1][1]))
        self.buffers_to_write = {}
        for renderer, _ in self.renderers.values():
            renderer.initializeGL(max_cells)
            write_funcs = renderer.set_buffers
            if write_funcs is not None:
                for field, func in write_funcs.items():
                    if field in self.buffers_to_write:
                        self.buffers_to_write[field] += [func]
                    else:
                        self.buffers_to_write[field] = [func]

    def resizeGL(self, w: int, h: int) -> None:
        """
        @brief from @c QtWidgets.QOpenGLWidget.
        """
        gl.glViewport(0, 0, w, h)
        self.projection.setToIdentity()
        self.projection.perspective(
            FOV, self.width() / self.height(), NEAR_DIST, FAR_DIST
        )
        self.projection_inv = self.projection.inverted()[0]

    def paintGL(self) -> None:
        """
        @brief from @c QtWidgets.QOpenGLWidget.
        """
        gl.glClearColor(0.5, 0.5, 0.5, 0.0)
        gl.glClear(int(gl.GL_COLOR_BUFFER_BIT) | int(gl.GL_DEPTH_BUFFER_BIT))
        gl.glEnable(gl.GL_DEPTH_TEST)

        # send cell data to renderers that need it
        for renderer, _ in self.renderers.values():
            if self.sim is not None:
                renderer.paintGL(
                    self.projection,
                    self.view,
                    self.projection_inv,
                    self.view_inv,
                    len(self.sim.cell_states),
                )
            else:
                renderer.paintGL(
                    self.projection, self.view, self.projection_inv, self.view_inv
                )

        if self.sim is not None:
            self.sim.write_to_buffer(self.buffers_to_write)

    ## @}

    ## @name PyQt5 Signals and slots
    # Interacts with PyGLGUI.ui components
    # @{

    ## @brief Signal for "Run" button to match current state.
    toggleRunSignal = pyqtSignal(bool)

    ## @brief Signal for "Save Pickles" button to match configured state.
    setSavePicklesToggle = pyqtSignal(bool)

    ## @brief Signal for selecting cells.
    selectedCell = pyqtSignal(str)  # emit selected cell info

    @pyqtSlot(bool)
    def toggleRun(self, run: bool) -> None:
        """
        @brief Handle for 'Run' button.
        """
        self.app_state["run"] = run
        if self.sim is not None:
            self.sim.paused = not run

    @pyqtSlot(bool)
    def toggleSavePickles(self, save: bool) -> None:
        """
        @brief Handle for 'Save pickles' button.
        """
        del save
        # if self.sim is not None:
        #     self.sim.setSaveOutput(save)

    @pyqtSlot()
    def reset(self) -> None:
        """
        @brief Handle for 'Reset' button.
        """
        sim = Simulator(
            self.modName,
            dt=DT,
            clPlatformNum=self.clPlatformNum,
            clDeviceNum=self.clDeviceNum,
            is_gui=True,
        )
        self._init_sim(sim)

        self.app_state["run"] = False
        self.toggleRunSignal.emit(False)
        self.update()

    @pyqtSlot()
    def loadGeometry(self) -> None:
        """
        @brief Handle for 'Load Geometry' button.
        @details Opens a file dialog chooser for a pickle file.
        @warning Pickles broken for now.
        """
        QMessageBox.warning(
            self,
            "Warning",
            "Pickles not yet supported on this branch.",
            QMessageBox.StandardButton.Ok,
        )

    @pyqtSlot()
    def loadPickle(self) -> None:
        """
        @brief Handle for 'Load Pickle' button.
        @details Opens a file dialog chooser for a pickle file.
        @warning Pickles broken for now.
        """
        QMessageBox.warning(
            self,
            "Warning",
            "Pickles not yet supported on this branch.",
            QMessageBox.StandardButton.Ok,
        )

    @pyqtSlot()
    def load(self) -> None:
        """
        @brief Handle for 'Load Model' button.
        @details Open a file dialog chooser for a python module.
        """
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        qs, _ = QFileDialog.getOpenFileName(
            self, "Load Python module", "", "*.py", options=options
        )
        if qs:
            modfile = str(qs)
            print(modfile)
            self.loadModelFile(modfile)

    def loadModelFile(self, modname: str) -> None:
        """
        @brief Bypass file dialog for choosing model.
        @param modname Absolute path to module python file.
        """
        if self._getOpenCLPlatDev():
            self.modName = modname
            try:
                sim = Simulator(
                    self.modName,
                    dt=DT,
                    clPlatformNum=self.clPlatformNum,
                    clDeviceNum=self.clDeviceNum,
                    is_gui=True,
                )
            except ModuleTypeError:
                QMessageBox.critical(
                    self,
                    "Error",
                    f"No UserModule found in {modname}.",
                    QMessageBox.StandardButton.Ok,
                )
                return

            self._init_sim(sim)
            self.update()

    ## @}

    ## @name Private methods
    # Manage simulation running.
    # @{

    def _init_sim(self, sim: Simulator) -> None:
        """
        @brief Attach configured simulator object.
        """
        if self.sim is not None:
            del self.sim
        self.sim = sim

        # Make GUI button match simulator state for saving pickles
        # self.setSavePicklesToggle.emit(sim.saveOutput)
        # print("saveOutput", sim.saveOutput)

        self.renderers = self.default_renderers
        for renderer in sim.user_module.renderers:
            self.renderers |= {renderer.__class__.__name__: (renderer, 0)}

        self.initialize_renderers()

    def _getOpenCLPlatDev(self) -> bool:
        """
        @brief Dialog for choosing compute platform.
        @details Opens an input chooser if multiple platforms detected. If only
                1 is available then default to that, skipping dialog.
        @return Indicates whether or not a platform was selected. Automatically
                updates internal platform id.
        """
        # Pop dialogs to get user to choose OpenCL platform
        platforms = cl.get_platforms()
        devices = platforms[self.clPlatformNum].get_devices()

        platlist = [str(p.name) for p in platforms]
        platdict = dict(list(zip(platlist, list(range(len(platlist))))))
        devlist = [str(d.name) for d in devices]
        devdict = dict(list(zip(devlist, list(range(len(devlist))))))

        if len(platlist) == 1:
            self.clPlatformNum = 0
            self.clDeviceNum = 0
            return True

        qsPlatformName, ok = QInputDialog.getItem(
            self,
            "Choose OpenCL platform",
            "Available platforms:",
            platlist,
            editable=False,
        )
        qsDeviceName, ok = QInputDialog.getItem(
            self, "Choose OpenCL device", "Available devices:", devlist, editable=False
        )

        if not ok:
            print("You didn't select a OpenCL platform...")
            return False
        else:
            self.clPlatformNum = platdict[qsPlatformName]
            self.clDeviceNum = devdict[qsDeviceName]
            return True

    class UpdateSim(QThread):
        """
        @brief QThread for managing simulation handling.
        @details The goal is to not block the render loop while waiting on sim
                steps. Framerate drops can still occur though...
        """

        def __init__(
            self,
            sim: Simulator,
            toggleRunSignal: pyqtBoundSignal,
            state: dict[str, bool],
            parent: Optional[QObject] = None,
        ) -> None:
            self.sim = sim
            self.toggleRunSignal = toggleRunSignal
            self.state = state
            super().__init__(parent)
            self.counter = 0

        def run(self) -> None:
            """
            @brief Runs the simulation step in a separate thread.
            @details Checks if sim has paused itself so that GUI mirrors it.
            """
            if self.sim.paused and self.state["run"]:
                self.toggleRunSignal.emit(False)
                return
            if self.state["run"]:
                self.sim.step(DT)

    def _animate(self) -> None:
        # triggers repaint
        self.update()

        # do regular updates during event loop
        dt = time.time() - self.last_time
        self._handle_keys(dt)
        self._get_fps(dt)
        self.last_time = time.time()

        # step sim only if at least DT real time has passed
        self.delta += dt
        if self.delta >= DT:
            if self.sim_thread is None or self.sim_thread.isFinished():
                if self.sim is not None:
                    self.sim_thread = self.UpdateSim(
                        self.sim, self.toggleRunSignal, self.app_state
                    )
                    self.sim_thread.start()
            self.delta -= DT

    ## @}


# Hide full stack if exception occurs in UpdateSim thread
def suppress_errors(exctype, value, traceback):
    if value is PermissionError:
        return
    sys.__excepthook__(exctype, value, traceback)
    exit(1)


sys.excepthook = suppress_errors
