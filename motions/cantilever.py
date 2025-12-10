from manim import *
import numpy as np


class CantileverBeam(VGroup):
    """A cantilever beam with fixed support at left end."""
    
    def __init__(self, length=6, height=0.4, **kwargs):
        super().__init__(**kwargs)
        self.length = length
        self.height = height
        
        # Create horizontal layers with alternating colors
        self.num_layers = 10
        self.layers = []
        layer_height = height / self.num_layers
        colors = [BLUE_D, BLUE_C]
        
        for i in range(self.num_layers):
            # Distance from neutral axis (center of beam)
            y_from_center = -height/2 + (i + 0.5) * layer_height
            
            layer = Rectangle(
                width=length,
                height=layer_height,
                color=colors[i % 2],
                fill_opacity=0.92,
                stroke_width=0
            )
            layer.move_to([0, y_from_center, 0])
            layer.y_from_center = y_from_center
            
            self.layers.append(layer)
            self.add(layer)
    
    def get_deflection_at_x(self, x, max_deflection):
        """Cubic deflection for cantilever with load at free end.
        
        Deflection shape: y(x) = -delta_max * (3*x^2/L^2 - 2*x^3/L^3)
        This matches the standard cantilever deflection curve.
        """
        L = self.length
        normalized_x = (x + L/2) / L  # normalize to [0, 1]
        normalized_x = max(0, min(1, normalized_x))  # clamp to [0, 1]
        return -max_deflection * (3 * normalized_x**2 - 2 * normalized_x**3)
    
    def get_slope_at_x(self, x, max_deflection):
        """Slope of deflection curve (derivative)."""
        L = self.length
        normalized_x = (x + L/2) / L
        normalized_x = max(0, min(1, normalized_x))
        # d/dx of -delta * (3*x^2/L^2 - 2*x^3/L^3) = -delta * (6*x/L^2 - 6*x^2/L^3)
        return -max_deflection * (6 * normalized_x / L - 6 * normalized_x**2 / L)
    
    def animate_deflection(self, max_deflection=0.5, run_time=2):
        """Animate proper beam bending using an internal alpha timeline."""
        return UpdateFromAlphaFunc(self, lambda m, a: self.apply_deflection(a * max_deflection), run_time=run_time)

    def apply_deflection(self, current_deflection: float, samples: int = 300):
        """Apply a given deflection value to all layers (no animation wrapper)."""
        layer_height = self.height / self.num_layers
        num_points = max(40, int(samples))

        for layer in self.layers:
            y_from_center = layer.y_from_center

            top_points = []
            bottom_points = []

            for i in range(num_points):
                t = i / (num_points - 1)
                x = -self.length / 2 + t * self.length

                # Neutral axis deflection and slope
                deflection = self.get_deflection_at_x(x, current_deflection)
                slope = self.get_slope_at_x(x, current_deflection)

                # Unit perpendicular to curve
                tangent_mag = np.sqrt(1 + slope**2)
                perp_x = -slope / tangent_mag
                perp_y = 1 / tangent_mag

                # Neutral axis position
                neutral_x = x
                neutral_y = deflection

                # This layer's center position
                center_x = neutral_x + perp_x * y_from_center
                center_y = neutral_y + perp_y * y_from_center

                # Top and bottom of layer
                half_h = layer_height / 2
                top_points.append([center_x + perp_x * half_h, center_y + perp_y * half_h, 0])
                bottom_points.append([center_x - perp_x * half_h, center_y - perp_y * half_h, 0])

            # Close the polygon
            all_points = top_points + bottom_points[::-1]
            layer.set_points_as_corners(all_points + [all_points[0]])

    def add_deflection_updater(self, tracker: ValueTracker, samples: int = 300):
        """Attach an updater that keeps the beam shape synced to `tracker`.

        Returns the updater function so it can be removed later if needed.
        """
        def _upd(mob):
            self.apply_deflection(tracker.get_value(), samples=samples)
        self.add_updater(_upd)
        # Ensure initial state is consistent with tracker value (no visual jump)
        self.apply_deflection(tracker.get_value(), samples=samples)
        return _upd
    
    def get_top_surface_point_at_x(self, x, current_deflection):
        """Get the y-coordinate of top surface at position x"""
        deflection = self.get_deflection_at_x(x, current_deflection)
        slope = self.get_slope_at_x(x, current_deflection)
        
        # Unit perpendicular
        tangent_mag = np.sqrt(1 + slope**2)
        perp_x = -slope / tangent_mag
        perp_y = 1 / tangent_mag
        
        # Top surface is neutral axis + height/2 in perpendicular direction
        top_x = x + perp_x * self.height / 2
        top_y = deflection + perp_y * self.height / 2
        
        return top_x, top_y

    def get_corner_position(self, x, current_deflection):
        """Get the corner position (top-left/right edge) of the beam at position x."""
        deflection = self.get_deflection_at_x(x, current_deflection)
        slope = self.get_slope_at_x(x, current_deflection)
        
        # Unit perpendicular
        tangent_mag = np.sqrt(1 + slope**2)
        perp_x = -slope / tangent_mag
        perp_y = 1 / tangent_mag
        
        # Corner is at neutral axis + height/2 in perpendicular direction
        corner_x = x + perp_x * self.height / 2
        corner_y = deflection + perp_y * self.height / 2
        
        return corner_x, corner_y


class CantileverLoadDeflection(Scene):
    def construct(self):
        # Create cantilever beam
        L = 6
        beam = CantileverBeam(length=L, height=0.4)
        # Ensure starting geometry uses the same path logic to avoid jumps
        beam.apply_deflection(0.0)
        self.play(FadeIn(beam), run_time=0.6)
        self.wait(0.5)

        # Point load at free end (x = L/2)
        delta_max = 0.8
        defl = ValueTracker(0.0)
        beam.add_deflection_updater(defl, samples=360)

        def free_end_top_pos():
            x, y = beam.get_top_surface_point_at_x(L / 2, defl.get_value())
            return np.array([x, y, 0.0])

        # Always keep the arrow in contact with the beam's top surface at free end
        arrow = always_redraw(
            lambda: Arrow(
                start=free_end_top_pos() + 1.4 * UP,
                end=free_end_top_pos(),
                color=RED,
                buff=0,
                stroke_width=7,
                tip_length=0.18,
            )
        )
        label = always_redraw(lambda: MathTex("P", color=RED, font_size=52).next_to(free_end_top_pos() + 1.4 * UP, UP, buff=0.15))

        # Fixed support at left end (x = -L/2)
        tri_width = 0.52
        tri_height = 0.18
        fixed_x = -L / 2
        fixed_y = -0.2  # fixed height
        
        # Create fixed support with triangular seat and ground hatching
        fixed_support = VGroup(
            # Triangle pointing up (fixed support)
            Polygon(
                [fixed_x - tri_width/2, fixed_y - tri_height, 0],
                [fixed_x + tri_width/2, fixed_y - tri_height, 0],
                [fixed_x, fixed_y, 0],  # apex
                color=WHITE,
                fill_opacity=1,
                stroke_width=2,
                stroke_color=WHITE
            ),
            # Base plate
            Rectangle(
                width=0.8, 
                height=0.12, 
                color=GRAY,
                fill_opacity=1
            ).move_to([fixed_x, fixed_y - 0.52, 0]),
            # Ground hatching
            VGroup(*[
                Line(
                    [fixed_x - 0.36 + i*0.11, fixed_y - 0.6, 0],
                    [fixed_x - 0.25 + i*0.11, fixed_y - 0.78, 0],
                    color=WHITE,
                    stroke_width=2
                ) for i in range(8)
            ])
        )

        self.add(fixed_support)

        self.play(Create(arrow), FadeIn(label), run_time=0.6)
        self.wait(0.2)

        # Animate the deflection smoothly via shared tracker
        self.play(defl.animate.set_value(delta_max), run_time=2.5, rate_func=smooth)
        # Freeze the current shape (optional): remove updaters to reduce CPU in long waits
        beam.clear_updaters()
        self.wait(1)
        
        # Formula for cantilever deflection
        formula = MathTex(
            r"\delta_{max} = \frac{PL^3}{3EI}",
            font_size=42,
            color=WHITE
        ).to_corner(UR).shift(DOWN * 0.3)
        
        self.play(Write(formula), run_time=1)
        self.wait(2)
