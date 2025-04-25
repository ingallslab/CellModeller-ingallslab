"""
@brief Class to handle initialization and dispatch draw calls.
"""

import importlib.resources
from abc import ABC, abstractmethod
from collections.abc import Callable

import numpy as np
import OpenGL.GL as gl
from PyQt5.QtGui import (
    QMatrix4x4,
    QOpenGLShader,
    QOpenGLShaderProgram,
    QOpenGLVertexArrayObject,
)

## @brief Type signature for buffer copy function.
WriteFunc = Callable[[memoryview], None]

## @brief Conversion between numpy format specifier and OpenGL types.
GL_TYPES = {
    "i": gl.GL_INT,
    "f": gl.GL_FLOAT,
}


class Renderer(ABC):
    """
    @brief Abstract base class for OpenGL renderers.
    @details Compiles shaders and automatically sets global uniforms. Provides a
            utility function to generate vertex attribute array from dtype.
    """

    ## @brief Optional cell attributes that the renderer requires.
    cell_attrs: np.dtype | None = None
    ## @brief Optional functions to copy a memoryview into a VBO.
    set_buffers: dict[str, WriteFunc] | None = None

    ## @brief OpenGL shader program.
    # @details Required to be compiled before binding attribute locations.
    shader_program: QOpenGLShaderProgram
    ## @brief Vertex array for this class.
    # @details Classes manage their own vertex buffers if needed. Assumed that
    #       each renderer only needs one VAO.
    vertex_array: QOpenGLVertexArrayObject

    ## @name Global uniform locations
    # Set by Renderer after shader is compiled.
    # @{

    _view_loc: int
    _view_inv_loc: int
    _proj_loc: int
    _proj_inv_loc: int

    # @}

    def compile_shader(self) -> None:
        """
        @brief Attach and configure shader before handing off to subclass.
        @details Will look for glsl programs in the same directory as the
                subclassing renderer and with the same name. File extensions
                are as follows:
                    - .vert for vertex shaders
                    - .geom for geometry shaders
                    - .frag for fragment shaders
        """
        self.shader_program = QOpenGLShaderProgram()

        name = self.__class__.__name__
        module = self.__class__.__module__

        vert = name + ".vert"
        try:
            with importlib.resources.open_text(module, vert) as vert_src:
                if not self.shader_program.addShaderFromSourceCode(
                    QOpenGLShader.ShaderTypeBit.Vertex, vert_src.read()
                ):
                    raise RuntimeError(f"Could not compile {vert} in {module}")
        except FileNotFoundError:
            pass
        geom = name + ".geom"
        try:
            with importlib.resources.open_text(module, geom) as geom_src:
                if not self.shader_program.addShaderFromSourceCode(
                    QOpenGLShader.ShaderTypeBit.Geometry, geom_src.read()
                ):
                    raise RuntimeError(f"Could not compile {geom} in {module}")
        except FileNotFoundError:
            pass
        frag = name + ".frag"
        try:
            with importlib.resources.open_text(module, frag) as frag_src:
                if not self.shader_program.addShaderFromSourceCode(
                    QOpenGLShader.ShaderTypeBit.Fragment, frag_src.read()
                ):
                    raise RuntimeError(f"Could not compile {frag} in {module}")
        except FileNotFoundError:
            pass

        if not self.shader_program.link():
            raise RuntimeError(f"Could not link shader program in {module}")

        self._view_loc = self.shader_program.uniformLocation("view")
        self._view_inv_loc = self.shader_program.uniformLocation("viewInv")
        self._proj_loc = self.shader_program.uniformLocation("proj")
        self._proj_inv_loc = self.shader_program.uniformLocation("projInv")

    def get_buffer_attrs(
        self, offset: int, max_cells: int
    ) -> tuple[int, dict[str, tuple[int, int, int, int]]]:
        """
        @brief Get VAO parameters from dtype.
        @param offset In bytes to start from.
        @param max_cells Required to match size of CellArrays array.
        @return Tuple of the new offset including the size of all added and
                dict mapping attribute names to tuples of glVertexAttribPointer
                arugments, that is: type, offset, tupleSize, stride. Attrib
                location determined in the subclass since attribs can only be
                set after the VBO has been allocated, which usually requires
                the calculated offset.
        """
        buffer_attributes = {}
        if self.cell_attrs is not None and self.cell_attrs.fields is not None:
            for field, (dtype, *_) in self.cell_attrs.fields.items():
                if dtype.subdtype:
                    subdtype, shape = dtype.subdtype
                else:
                    subdtype, shape = (dtype, (1,))
                buffer_attributes[field] = (
                    GL_TYPES[subdtype.char],
                    offset,
                    subdtype.itemsize,
                    subdtype.itemsize * shape[0],
                )
                offset += dtype.itemsize * max_cells
        return (offset, buffer_attributes)

    def initializeGL(self, max_cells: int) -> None:
        """
        @brief Called by PyGLCMViewer once OpenGL context is ready.
        @details Compiles shader and then calls init_renderer().
        """
        self.compile_shader()
        self.init_renderer(max_cells)

    @abstractmethod
    def init_renderer(self, max_cells: int) -> None:
        """
        @brief Called once the shader has been compiled.
        @details Must be implemented in subclasses and initializes additional
                renderer parameters, creates VBOs, etc.
        @param max_cells Passed on to get_buffer_attrs().
        """
        ...

    def paintGL(
        self,
        view: QMatrix4x4,
        proj: QMatrix4x4,
        view_inv: QMatrix4x4,
        proj_inv: QMatrix4x4,
        vertex_count: int = 0,
    ) -> None:
        """
        @brief Called by PyGLCMViewer once every render cycle.
        @details Bind and unbind shader program and VAO, setting global uniforms
                and calling the subclass draw() function in between.
        """
        self.shader_program.bind()
        self.vertex_array.bind()

        self.shader_program.setUniformValue(self._view_loc, view)
        self.shader_program.setUniformValue(self._proj_loc, proj)
        self.shader_program.setUniformValue(self._view_inv_loc, view_inv)
        self.shader_program.setUniformValue(self._proj_inv_loc, proj_inv)

        self.draw(vertex_count)

        self.vertex_array.release()
        self.shader_program.release()

    @abstractmethod
    def draw(self, vertex_count: int) -> None:
        """
        @brief Main rendering logic goes here.
        @param vertex_count Equal to last index in cell array.
        """
        ...
