"""
@brief Class for background grid in GUI.
"""

import OpenGL.GL as gl
from PyQt5.QtGui import QVector3D

from CellModeller.GUI.Renderers.Renderer import Renderer


class GridRenderer(Renderer):
    """
    @brief Subclasses Renderer to handle OpenGL boilerplate.
    @details Parameters that can be configured within the class:
            - minTick the spacing of grid lines
            - majTick the spacing of bold grid lines
            - linecolor of off-axis lines
    """

    def __init__(self) -> None:
        self.uniforms: dict[str, int | QVector3D] = {
            "minTick": 5,
            "majTick": 20,
            "linecolor": QVector3D(0.8, 0.8, 0.8),
        }

    def init_renderer(self, max_cells: int) -> None:
        del max_cells

    def draw(self, vertex_count: int) -> None:
        del vertex_count
        gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
        gl.glEnable(gl.GL_BLEND)
        gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)
        gl.glDisable(gl.GL_CULL_FACE)

        for name, val in self.uniforms.items():
            loc = self.shader_program.uniformLocation(name)
            self.shader_program.setUniformValue(loc, val)

        # vertices hard-coded in vertex shader
        gl.glDrawArrays(gl.GL_TRIANGLES, 0, 6)


"""
Adapted from Dube, M. (2020, January 5). How to make an infinite grid.
https://asliceofrendering.com/scene%20helper/2020/01/05/InfiniteGrid/.
"""
