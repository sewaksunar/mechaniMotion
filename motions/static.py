from manim import *
import numpy as np

class MySquare(Square):
    @override_animation(FadeIn)
    def _fade_in_override(self, **kwargs):
        return Create(self, **kwargs)

class OverrideAnimationExample(Scene):
    def construct(self):
        self.play(FadeIn(MySquare()))

class PinSupportBar(VGroup):
    def __init__(self, length=4, height=0.5, **kwargs):
        super().__init__(**kwargs)
        self.length = length
        self.height = height
        self.bar = Rectangle(width=length, height=height, color=BLUE)
        pin1 = Triangle(color=RED).scale(0.3).next_to(self.bar.get_left(), DOWN, buff=0)
        pin2 = Triangle(color=RED).scale(0.3).next_to(self.bar.get_right(), DOWN, buff=0)
        pin2.shift(DOWN * 0.2)
        pin1.shift(DOWN * 0.2)
        self.add(self.bar, pin1, pin2)

    def deflected_curve(self, max_deflection=0.3, samples=60, color=YELLOW):
        """
        Returns a smooth curve representing the deflected shape of a simply
        supported bar under a point load at the center. We approximate the
        deflection shape with y(x) = delta_max * sin(pi * x / L) for visual clarity.

        max_deflection: peak deflection at midspan in scene units
        samples: number of points used to draw the curve
        """
        L = self.length
        left = self.bar.get_left()
        # Build local coordinate system along the bar width
        pts = []
        for i in range(samples + 1):
            t = i / samples
            x = t * L
            # sine mode approximation (matches zero deflection at supports)
            y = -max_deflection * np.sin(np.pi * x / L)
            # position relative to left support
            pts.append(left + RIGHT * x + UP * y)
        curve = VMobject(color=color, stroke_width=6)
        curve.set_points_smoothly(pts)
        return curve

    def animate_deflection(self, max_deflection=0.3, run_time=2):
        curve = self.deflected_curve(max_deflection=max_deflection)
        return Transform(self.bar, curve, run_time=run_time)
class PinSupportBarExample(Scene):

    def construct(self):
        pin_support_bar = PinSupportBar()
        self.play(FadeIn(pin_support_bar))
        self.wait(1)
        transform_bar = PinSupportBar(length=6, height=0.5)
        self.play(Transform(pin_support_bar, transform_bar))

from manim import *

class PinSupportBar(VGroup):
    def __init__(self, length=6, height=0.4, **kwargs):
        super().__init__(**kwargs)
        self.length = length
        self.height = height
        
        # Store strip information
        self.num_strips = 20
        self.strips = []
        strip_height = height / self.num_strips
        
        # Create horizontal strips with alternating colors
        colors = [BLUE_D, BLUE_C]
        
        for i in range(self.num_strips):
            # Y position of strip center (from bottom to top)
            y_center = -height/2 + (i + 0.5) * strip_height
            
            # Create strip as rectangle
            strip = Rectangle(
                width=length,
                height=strip_height,
                color=colors[i % 2],
                fill_opacity=0.9,
                stroke_width=0
            )
            strip.move_to([0, y_center, 0])
            strip.y_offset = y_center  # Store original y position
            
            self.strips.append(strip)
            self.add(strip)
        
        # Left support - fixed pin
        left_x = -length/2
        self.left_support = VGroup(
            Polygon(
                [left_x - 0.3, -height/2 - 0.4, 0],
                [left_x + 0.3, -height/2 - 0.4, 0],
                [left_x, -height/2, 0],
                color=WHITE,
                fill_opacity=1,
                stroke_width=2
            ),
            Rectangle(
                width=0.8, 
                height=0.15, 
                color=GRAY,
                fill_opacity=1
            ).move_to([left_x, -height/2 - 0.55, 0]),
            VGroup(*[
                Line(
                    [left_x - 0.35 + i*0.12, -height/2 - 0.65, 0],
                    [left_x - 0.23 + i*0.12, -height/2 - 0.85, 0],
                    color=WHITE,
                    stroke_width=2
                ) for i in range(7)
            ])
        )
        
        # Right support - roller
        right_x = length/2
        self.right_support = VGroup(
            Circle(
                radius=0.2,
                color=WHITE,
                fill_opacity=1,
                stroke_width=2
            ).move_to([right_x, -height/2 - 0.28, 0]),
            Rectangle(
                width=0.8, 
                height=0.15, 
                color=GRAY,
                fill_opacity=1
            ).move_to([right_x, -height/2 - 0.55, 0]),
            VGroup(
                Circle(radius=0.09, color=GRAY_BROWN, fill_opacity=1).move_to([right_x - 0.22, -height/2 - 0.65, 0]),
                Circle(radius=0.09, color=GRAY_BROWN, fill_opacity=1).move_to([right_x + 0.22, -height/2 - 0.65, 0])
            ),
            VGroup(*[
                Line(
                    [right_x - 0.35 + i*0.12, -height/2 - 0.7, 0],
                    [right_x - 0.23 + i*0.12, -height/2 - 0.9, 0],
                    color=WHITE,
                    stroke_width=2
                ) for i in range(7)
            ])
        )
        
        self.add(self.left_support, self.right_support)
    
    def get_deflection_at_x(self, x, max_deflection):
        """Calculate deflection using parabolic curve"""
        normalized_x = 2 * x / self.length
        return -max_deflection * (1 - normalized_x**2)
    
    def animate_deflection(self, max_deflection=0.5, run_time=2):
        """Animate beam deflection with horizontal strips maintaining constant height"""
        strip_height = self.height / self.num_strips
        
        def apply_deflection(mob, alpha):
            current_deflection = alpha * max_deflection
            num_segments = 100  # Segments per strip for smooth curve
            
            for strip in self.strips:
                y_offset = strip.y_offset
                
                # Create curved path for this strip
                top_edge = []
                bottom_edge = []
                
                for j in range(num_segments + 1):
                    t = j / num_segments
                    x = -self.length/2 + t * self.length
                    
                    # Deflection at this x position
                    deflection = self.get_deflection_at_x(x, current_deflection)
                    
                    # Position of strip center at this x
                    y_center = y_offset + deflection
                    
                    # Top and bottom edges (maintain constant strip height)
                    top_edge.append([x, y_center + strip_height/2, 0])
                    bottom_edge.append([x, y_center - strip_height/2, 0])
                
                # Create closed polygon: top edge + reversed bottom edge
                points = top_edge + bottom_edge[::-1]
                strip.set_points_as_corners(points + [points[0]])
        
        return UpdateFromAlphaFunc(self, apply_deflection, run_time=run_time)


class CenterLoadDeflection(Scene):
    def construct(self):
        # Create beam
        L = 6
        bar = PinSupportBar(length=L, height=0.4)
        self.play(FadeIn(bar), run_time=1)
        self.wait(0.5)

        # Point load at center
        arrow = Arrow(
            start=[0, 1.6, 0], 
            end=[0, 0.25, 0], 
            color=RED, 
            buff=0, 
            stroke_width=7,
            max_tip_length_to_length_ratio=0.15
        )
        label = MathTex("P", color=RED, font_size=52).next_to(arrow.get_start(), UP, buff=0.15)
        
        self.play(GrowArrow(arrow), FadeIn(label), run_time=0.8)
        self.wait(0.4)

        # Animate deflection
        delta_max = 0.8
        
        def update_arrow(mob, alpha):
            current_deflection = alpha * delta_max
            y_top = 0.2 + current_deflection
            mob.put_start_and_end_on([0, y_top + 1.4, 0], [0, y_top, 0])
        
        label.add_updater(lambda m: m.next_to(arrow.get_start(), UP, buff=0.15))
        
        self.play(
            bar.animate_deflection(max_deflection=delta_max, run_time=2.5),
            UpdateFromAlphaFunc(arrow, update_arrow, run_time=2.5),
            rate_func=smooth
        )
        
        label.clear_updaters()
        self.wait(1)
        
        # Formula
        formula = MathTex(
            r"\delta_{max} = \frac{PL^3}{48EI}",
            font_size=42,
            color=WHITE
        ).to_corner(UR).shift(DOWN * 0.3)
        
        self.play(Write(formula), run_time=1)
        self.wait(2)