from manim import *


def main():
    print("Hello from manimations!")


if __name__ == "__main__":
    main()


class CreateCircle(Scene):
    def construct(self):
        circle = Circle()  # create a circle
        circle.set_fill(PINK, opacity=0.5)  # set the color and transparency
        self.play(Create(circle))  # show the circle on screen