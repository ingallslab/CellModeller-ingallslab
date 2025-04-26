"""
@brief Class for rendering default rod-shaped cells in application.
"""

import numpy as np
import OpenGL.GL as gl
from PyQt5.QtGui import QOpenGLBuffer

from CellModeller.GUI.Renderers.Renderer import Renderer, WriteFunc


class BacillusRenderer(Renderer):
    """
    @brief Renderer handles OpenGL boilerplate except VBO management.
    @details Cell attributes defined by this class can be used by all other
            modules, these include:
                - radius
                - length
                - pos
                - dir
                - color
    @param wireframe Whether or not to render in wireframe mode. (it looks cool)
    """

    def __init__(self, wireframe: bool = False) -> None:
        self.wireframe = wireframe
        self.cell_attrs = np.dtype(
            [
                ("radius", "f"),
                ("length", "f"),
                ("pos", "3f"),
                ("dir", "3f"),
                ("color", "3f"),
            ]
        )

    def write_from(self, offset: int) -> WriteFunc:
        """
        @brief Used to construct buffer copy functions in set_buffers.
        @details Necessarily bind/unbind the VBO outside of copying the buffer.
        """

        def bind_buffers(data: memoryview):
            self.vertex_buffer.bind()
            self.vertex_buffer.write(offset, data, data.nbytes)  # type: ignore[reportArgumentType]
            self.vertex_buffer.release()

        return bind_buffers

    def init_renderer(self, max_cells: int) -> None:
        self.block_start = max_cells
        size, buffer_attributes = self.get_buffer_attrs(self.block_start, max_cells)
        self.set_buffers = {
            "BacillusRenderer": self.write_from(self.block_start),
            "_active_flag": self.write_from(0),
        }

        self.vertex_buffer = QOpenGLBuffer(QOpenGLBuffer.Type.VertexBuffer)
        self.vertex_buffer.create()

        self.vertex_buffer.bind()
        self.vertex_buffer.setUsagePattern(QOpenGLBuffer.UsagePattern.DynamicDraw)
        self.vertex_buffer.allocate(size)

        # required so we don't render deleted cells...
        loc = self.shader_program.attributeLocation("activeFlag")
        self.shader_program.enableAttributeArray(loc)
        self.shader_program.setAttributeBuffer(loc, gl.GL_UNSIGNED_BYTE, 0, 1, 1)

        for name, attributes in buffer_attributes.items():
            loc = self.shader_program.attributeLocation(name)
            self.shader_program.enableAttributeArray(loc)
            self.shader_program.setAttributeBuffer(loc, *attributes)

        self.vertex_buffer.release()

        self.outline_loc = self.shader_program.uniformLocation("outline")

    def draw(self, vertex_count: int) -> None:
        if self.wireframe:
            gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_LINE)
        else:
            gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)
        gl.glFrontFace(gl.GL_CCW)
        gl.glEnable(gl.GL_CULL_FACE)

        # very neat trick i used here cull the opposite face to render the
        # "inside" faces of the outline, so that it doesn't obscure the capsule!
        gl.glCullFace(gl.GL_FRONT)
        self.shader_program.setUniformValue(self.outline_loc, True)
        gl.glDrawArrays(gl.GL_POINTS, 0, vertex_count)

        gl.glCullFace(gl.GL_BACK)
        self.shader_program.setUniformValue(self.outline_loc, False)
        gl.glDrawArrays(gl.GL_POINTS, 0, vertex_count)
