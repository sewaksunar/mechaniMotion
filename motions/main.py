def main():
    print("Hello from motions!")


if __name__ == "__main__":
    main()

from manim import *



from manim import *
import numpy as np

class StressTransformation(Scene):
    def construct(self):
        # ============ CONFIGURABLE PARAMETERS ============
        # Rotation angle
        ROTATION_ANGLE = PI/4  # Change this: PI/4 for 45°, PI/3 for 60°, etc.
        
        # Element properties
        ELEMENT_SIZE = 1.5
        ELEMENT_COLOR = BLUE
        TRANSFORMED_COLOR = PURPLE
        
        # Arrow properties
        ARROW_LENGTH = 0.6
        NORMAL_STRESS_COLOR = RED
        SHEAR_STRESS_COLOR = YELLOW
        
        # Axis properties
        ORIGINAL_AXIS_COLOR = GREEN
        ROTATED_AXIS_COLOR = ORANGE
        AXIS_LENGTH = 4
        
        # Animation timing
        ANIMATION_SPEED = 2.5
        WAIT_TIME = 1
        
        # Canvas layout
        VERTICAL_SHIFT = 0.5
        
        # ============ SETUP ============
        # Title
        title = Text("Stress Transformation Due to Rotation", font_size=32)
        title.to_edge(UP, buff=0.3)
        self.play(Write(title))
        self.wait(WAIT_TIME * 0.5)
        
        # Calculate angle in degrees for display
        angle_degrees = int(np.degrees(ROTATION_ANGLE))
        
        # ============ COORDINATE SYSTEM ============
        axes_original = Axes(
            x_range=[-2, 2, 1],
            y_range=[-2, 2, 1],
            x_length=AXIS_LENGTH,
            y_length=AXIS_LENGTH,
            axis_config={"color": ORIGINAL_AXIS_COLOR, "stroke_width": 2},
        )
        axes_original.shift(DOWN * VERTICAL_SHIFT)
        center_point = axes_original.get_center()
        
        # Original axis labels
        x_label = MathTex("x", color=ORIGINAL_AXIS_COLOR).scale(0.8)
        x_label.next_to(axes_original.x_axis.get_end(), RIGHT, buff=0.15)
        
        y_label = MathTex("y", color=ORIGINAL_AXIS_COLOR).scale(0.8)
        y_label.next_to(axes_original.y_axis.get_end(), UP, buff=0.15)
        
        # Show axes
        self.play(
            Create(axes_original),
            Write(VGroup(x_label, y_label))
        )
        self.wait(WAIT_TIME * 0.3)
        
        # ============ STRESS ELEMENT ============
        square = Square(side_length=ELEMENT_SIZE)
        square.set_fill(ELEMENT_COLOR, opacity=0.2)
        square.set_stroke(ELEMENT_COLOR, width=3)
        square.move_to(center_point)
        
        self.play(FadeIn(square))
        self.wait(WAIT_TIME * 0.3)
        
        # ============ STRESS VECTORS ============
        def create_stress_arrows(element, arrow_length, normal_color, shear_color):
            """Create stress arrows for a given element"""
            arrows = VGroup()
            
            # Normal stresses - σx (horizontal)
            arrows.add(Arrow(
                start=element.get_right() + LEFT * 0.1,
                end=element.get_right() + RIGHT * arrow_length,
                color=normal_color, buff=0, stroke_width=4,
                max_tip_length_to_length_ratio=0.25
            ))
            arrows.add(Arrow(
                start=element.get_left() + RIGHT * 0.1,
                end=element.get_left() + LEFT * arrow_length,
                color=normal_color, buff=0, stroke_width=4,
                max_tip_length_to_length_ratio=0.25
            ))
            
            # Normal stresses - σy (vertical)
            arrows.add(Arrow(
                start=element.get_top() + DOWN * 0.1,
                end=element.get_top() + UP * arrow_length,
                color=normal_color, buff=0, stroke_width=4,
                max_tip_length_to_length_ratio=0.25
            ))
            arrows.add(Arrow(
                start=element.get_bottom() + UP * 0.1,
                end=element.get_bottom() + DOWN * arrow_length,
                color=normal_color, buff=0, stroke_width=4,
                max_tip_length_to_length_ratio=0.25
            ))
            
            # Shear stresses - τxy
            arrows.add(Arrow(
                start=element.get_corner(UR) + LEFT * 0.1,
                end=element.get_corner(UR) + UP * arrow_length * 0.7,
                color=shear_color, buff=0, stroke_width=3,
                max_tip_length_to_length_ratio=0.3
            ))
            arrows.add(Arrow(
                start=element.get_corner(UR) + DOWN * 0.1,
                end=element.get_corner(UR) + RIGHT * arrow_length * 0.7,
                color=shear_color, buff=0, stroke_width=3,
                max_tip_length_to_length_ratio=0.3
            ))
            
            return arrows
        
        stress_arrows = create_stress_arrows(
            square, ARROW_LENGTH, NORMAL_STRESS_COLOR, SHEAR_STRESS_COLOR
        )
        
        # Stress labels
        sigma_x_label = MathTex(r"\sigma_x", color=NORMAL_STRESS_COLOR).scale(0.7)
        sigma_x_label.next_to(stress_arrows[0], RIGHT, buff=0.1)
        
        sigma_y_label = MathTex(r"\sigma_y", color=NORMAL_STRESS_COLOR).scale(0.7)
        sigma_y_label.next_to(stress_arrows[2], UP, buff=0.1)
        
        tau_label = MathTex(r"\tau_{xy}", color=SHEAR_STRESS_COLOR).scale(0.6)
        tau_label.next_to(stress_arrows[4], UR, buff=0.1)
        
        stress_labels = VGroup(sigma_x_label, sigma_y_label, tau_label)
        
        # Animate stress arrows
        self.play(
            *[GrowArrow(arrow) for arrow in stress_arrows],
            run_time=1.5
        )
        self.play(Write(stress_labels))
        self.wait(WAIT_TIME)
        
        # ============ ROTATION ANGLE ARC ============
        angle_arc = Arc(
            radius=1.2,
            start_angle=0,
            angle=ROTATION_ANGLE,
            color=YELLOW,
            stroke_width=3
        )
        angle_arc.move_arc_center_to(center_point)
        
        angle_label = MathTex(rf"\theta = {angle_degrees}°", color=YELLOW).scale(0.75)
        angle_label.next_to(angle_arc, RIGHT, buff=0.2)
        angle_label.shift(UP * 0.2)
        
        self.play(Create(angle_arc))
        self.play(Write(angle_label))
        self.wait(WAIT_TIME)
        
        # ============ ROTATED COORDINATE SYSTEM ============
        axes_rotated = Axes(
            x_range=[-2, 2, 1],
            y_range=[-2, 2, 1],
            x_length=AXIS_LENGTH,
            y_length=AXIS_LENGTH,
            axis_config={"color": ROTATED_AXIS_COLOR, "stroke_width": 2},
        )
        axes_rotated.move_to(center_point)
        axes_rotated.rotate(ROTATION_ANGLE, about_point=center_point)
        
        # Rotated axis labels
        x_prime_label = MathTex("x'", color=ROTATED_AXIS_COLOR).scale(0.8)
        y_prime_label = MathTex("y'", color=ROTATED_AXIS_COLOR).scale(0.8)
        
        # Group rotatable objects
        rotatable_group = VGroup(square, stress_arrows, stress_labels)
        
        # ============ ROTATION ANIMATION ============
        self.play(
            Create(axes_rotated),
            Rotate(rotatable_group, ROTATION_ANGLE, about_point=center_point),
            run_time=ANIMATION_SPEED
        )
        
        # Position rotated axis labels
        x_prime_label.next_to(axes_rotated.x_axis.get_end(), RIGHT, buff=0.15)
        y_prime_label.next_to(axes_rotated.y_axis.get_end(), UP, buff=0.15)
        
        self.play(Write(VGroup(x_prime_label, y_prime_label)))
        self.wait(WAIT_TIME)
        
        # ============ TRANSFORMED STRESS LABELS ============
        sigma_x_prime_label = MathTex(r"\sigma_{x'}", color=TRANSFORMED_COLOR).scale(0.7)
        sigma_x_prime_label.move_to(sigma_x_label)
        
        sigma_y_prime_label = MathTex(r"\sigma_{y'}", color=TRANSFORMED_COLOR).scale(0.7)
        sigma_y_prime_label.move_to(sigma_y_label)
        
        tau_prime_label = MathTex(r"\tau_{x'y'}", color=TRANSFORMED_COLOR).scale(0.6)
        tau_prime_label.move_to(tau_label)
        
        self.play(
            Transform(sigma_x_label, sigma_x_prime_label),
            Transform(sigma_y_label, sigma_y_prime_label),
            Transform(tau_label, tau_prime_label),
            square.animate.set_stroke(TRANSFORMED_COLOR, width=3),
            square.animate.set_fill(TRANSFORMED_COLOR, opacity=0.2),
        )
        self.wait(WAIT_TIME)
        
        # ============ TRANSFORMATION FORMULAS ============
        formula_box = Rectangle(
            width=12,
            height=2,
            color=WHITE,
            stroke_width=2
        )
        formula_box.to_edge(DOWN, buff=0.2)
        
        formulas = VGroup(
            MathTex(
                r"\sigma_{x'} = \frac{\sigma_x + \sigma_y}{2} + \frac{\sigma_x - \sigma_y}{2}\cos(2\theta) + \tau_{xy}\sin(2\theta)",
                font_size=28
            ),
            MathTex(
                r"\tau_{x'y'} = -\frac{\sigma_x - \sigma_y}{2}\sin(2\theta) + \tau_{xy}\cos(2\theta)",
                font_size=28
            )
        )
        formulas.arrange(DOWN, buff=0.3)
        formulas.move_to(formula_box.get_center())
        
        self.play(Create(formula_box))
        self.play(Write(formulas), run_time=2)
        
        self.wait(WAIT_TIME * 3)


class Shapes(Scene):
    def construct(self):
        circle = Circle()
        square = Square()
        triangle = Triangle()

        circle.shift(LEFT)
        square.shift(UP)
        triangle.shift(RIGHT)

        self.add(circle, square, triangle)

