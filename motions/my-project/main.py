from manim import *
class AxesTemplate(Scene):
    def construct(self):
        graph = Axes(
            x_range=[-1,10,1],
            y_range=[-1,10,1],
            x_length=9,
            y_length=6,
            axis_config={"include_tip":False}
        )
        labels = graph.get_axis_labels()
        self.add(graph, labels)
        circle = Circle()  # create a circle
        circle.set_fill(PINK, opacity=0.5)  # set the color and transparency
        self.add(circle)  # add the circle to the scene
