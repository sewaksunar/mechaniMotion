from manim import *
import numpy as np

class FourBarCurvature(MovingCameraScene):
    def construct(self):
        # Scaling: 1 mm -> SCALE scene units
        SCALE = 0.02
        VEL_SCALE = 0.004   # scales velocity arrows
        ACC_SCALE = 0.00025 # scales acceleration arrows

        # Geometry in millimeters (matches the textbook dimensions)
        RAO2 = 100.0                          # |O2A|
        RBA  = np.sqrt(2.0) * 100.0          # |AB|
        RCO4 = 100.0                          # |O4C|
        O2   = np.array([0.0, 0.0])
        O4   = np.array([-100.0, 0.0])

        # Reference posture (same as the figure)
        theta2_0 = 0.0
        theta3_0 = -3 * np.pi / 4  # -135°
        theta4_0 = -np.pi / 2      # -90°

        # First-order kinematic coefficients dθ3/dθ2 and dθ4/dθ2
        theta3_p = 0.5
        theta4_p = 0.5

        # Constant angular velocity of input link 2 (clockwise)
        omega2 = -10.0  # rad/s

        # --- Kinematics helpers ------------------------------------------------
        def geometry(theta2):
            """Return positions of A, B, C (in mm) and current θ3, θ4."""
            dtheta = theta2 - theta2_0
            theta3 = theta3_0 + theta3_p * dtheta
            theta4 = theta4_0 + theta4_p * dtheta

            A = O2 + RAO2 * np.array([np.cos(theta2), np.sin(theta2)])
            B = A  + RBA  * np.array([np.cos(theta3), np.sin(theta3)])
            C = O4 + RCO4 * np.array([np.cos(theta4), np.sin(theta4)])
            return A, B, C, theta3, theta4

        def kinematic_state(theta2):
            """Full kinematic state at a given input angle θ2."""
            A, B, C, theta3, theta4 = geometry(theta2)

            # Angular velocities of links 3 and 4 from first-order coefficients
            omega3 = theta3_p * omega2
            omega4 = theta4_p * omega2

            # Velocities (mm/s)
            vA = omega2 * RAO2 * np.array([-np.sin(theta2), np.cos(theta2)])
            vB = vA + omega3 * RBA * np.array([-np.sin(theta3), np.cos(theta3)])
            vC = omega4 * RCO4 * np.array([-np.sin(theta4), np.cos(theta4)])

            # Accelerations (mm/s^2) – input angular acceleration is zero here
            aA = (omega2 ** 2) * RAO2 * np.array([-np.cos(theta2), -np.sin(theta2)])
            aB = aA + (omega3 ** 2) * RBA * np.array([-np.cos(theta3), -np.sin(theta3)])
            aC = (omega4 ** 2) * RCO4 * np.array([-np.cos(theta4), -np.sin(theta4)])

            # Instantaneous center of curvature and radius for point B
            v = vB
            a = aB
            speed = np.linalg.norm(v)
            cross = v[0] * a[1] - v[1] * a[0]
            rho = speed ** 3 / abs(cross)  # radius of curvature (mm)

            # Unit normal pointing toward the center of curvature
            n = np.array([-v[1], v[0]])
            n /= np.linalg.norm(n)
            if np.dot(a, n) < 0:
                n = -n

            center = B + rho * n  # center of curvature position (mm)

            return dict(
                A=A, B=B, C=C,
                vA=vA, vB=vB, vC=vC,
                aA=aA, aB=aB, aC=aC,
                center=center, rho=rho,
            )

        def to_point(p_mm):
            """Convert 2D mm coordinates to Manim 3D point."""
            return np.array([p_mm[0], p_mm[1], 0.0]) * SCALE

        def vector_arrow(origin_mm, vec_mm, scale, color):
            """Arrow starting at origin_mm, pointing along vec_mm."""
            return Arrow(
                to_point(origin_mm),
                to_point(origin_mm + vec_mm * scale),
                max_tip_length_to_length_ratio=0.15,
                buff=0,
                stroke_width=3,
                color=color,
            )

        # Input angle tracker
        theta_tracker = ValueTracker(0.0)

        def state():
            return kinematic_state(theta_tracker.get_value())

        # --- Ground and fixed features ----------------------------------------
        ground = Line(
            to_point(np.array([-200.0, -100.0])),
            to_point(np.array([200.0, -100.0])),
            stroke_width=4,
            color=GREY_B,
        )

        O2_dot = Dot(to_point(O2), color=WHITE)
        O4_dot = Dot(to_point(O4), color=WHITE)

        O2_label = MathTex("O_2", font_size=28).next_to(O2_dot, DOWN + RIGHT, buff=0.1)
        O4_label = MathTex("O_4", font_size=28).next_to(O4_dot, DOWN + LEFT, buff=0.1)

        # --- Moving joints -----------------------------------------------------
        A_dot = always_redraw(lambda: Dot(to_point(state()["A"]), color=WHITE))
        B_dot = always_redraw(lambda: Dot(to_point(state()["B"]), color=YELLOW))
        C_dot = always_redraw(lambda: Dot(to_point(state()["C"]), color=WHITE))

        A_label = always_redraw(
            lambda: MathTex("A", font_size=28).next_to(A_dot, UP + RIGHT, buff=0.1)
        )
        B_label = always_redraw(
            lambda: MathTex("B", font_size=28).next_to(B_dot, DOWN + RIGHT, buff=0.1)
        )
        C_label = always_redraw(
            lambda: MathTex("C", font_size=28).next_to(C_dot, DOWN + LEFT, buff=0.1)
        )

        # --- Links -------------------------------------------------------------
        link2 = always_redraw(
            lambda: Line(
                O2_dot.get_center(),
                A_dot.get_center(),
                stroke_width=6,
                color=BLUE,
            )
        )
        link3 = always_redraw(
            lambda: Line(
                A_dot.get_center(),
                B_dot.get_center(),
                stroke_width=6,
                color=GREEN,
            )
        )
        link4 = always_redraw(
            lambda: Line(
                O4_dot.get_center(),
                C_dot.get_center(),
                stroke_width=6,
                color=PURPLE,
            )
        )
        # Horizontal bar between C and B (visual only, to mimic the book sketch)
        base_bar = always_redraw(
            lambda: Line(
                C_dot.get_center(),
                B_dot.get_center(),
                stroke_width=5,
                color=GREY_B,
            )
        )

        # Link labels
        link2_label = always_redraw(
            lambda: MathTex("2", font_size=26).move_to(
                (O2_dot.get_center() + A_dot.get_center()) / 2 + 0.25 * UP
            )
        )
        link3_label = always_redraw(
            lambda: MathTex("3", font_size=26).move_to(
                (A_dot.get_center() + B_dot.get_center()) / 2 + 0.25 * UP
            )
        )
        link4_label = always_redraw(
            lambda: MathTex("4", font_size=26).move_to(
                (O4_dot.get_center() + C_dot.get_center()) / 2 + 0.25 * LEFT
            )
        )

        # --- Path of point B ---------------------------------------------------
        path_B = TracedPath(
            lambda: B_dot.get_center(),
            stroke_color=YELLOW,
            stroke_width=3,
        )

        # --- Instantaneous center of curvature & curvature circle -------------
        curvature_center = always_redraw(
            lambda: Dot(to_point(state()["center"]), color=RED, radius=0.06)
        )
        curvature_circle = always_redraw(
            lambda: Circle(
                radius=state()["rho"] * SCALE,
                color=RED,
                stroke_width=2,
            ).move_to(curvature_center.get_center())
        )

        center_label = always_redraw(
            lambda: MathTex("C_{\\kappa}", font_size=24, color=RED).next_to(
                curvature_center, UP + LEFT, buff=0.1
            )
        )

        # --- Velocity and acceleration arrows at B ----------------------------
        vel_arrow = always_redraw(
            lambda: vector_arrow(
                state()["B"],
                state()["vB"],
                VEL_SCALE,
                color=BLUE,
            )
        )
        acc_arrow = always_redraw(
            lambda: vector_arrow(
                state()["B"],
                state()["aB"],
                ACC_SCALE,
                color=GREEN,
            )
        )

        vel_label = always_redraw(
            lambda: MathTex("\\vec{v}_B", font_size=24, color=BLUE).next_to(
                vel_arrow.get_end(), UP + RIGHT, buff=0.1
            )
        )
        acc_label = always_redraw(
            lambda: MathTex("\\vec{a}_B", font_size=24, color=GREEN).next_to(
                acc_arrow.get_end(), DOWN + RIGHT, buff=0.1
            )
        )

        # Angular-velocity annotation
        omega_text = MathTex(
            "\\omega_2 = 10\\,\\text{rad/s}~\\text{(cw)}",
            font_size=28,
        ).to_corner(UR)

        # Camera framing
        self.camera.frame.set_width(14)

        # Add everything to the scene
        self.add(
            ground,
            O2_dot, O4_dot, O2_label, O4_label,
            link2, link3, link4, base_bar,
            A_dot, B_dot, C_dot,
            A_label, B_label, C_label,
            link2_label, link3_label, link4_label,
            path_B,
            curvature_circle, curvature_center, center_label,
            vel_arrow, acc_arrow, vel_label, acc_label,
            omega_text,
        )

        self.wait(1)

        # Animate through a range of input angles
        theta_min = -0.7
        theta_max = 0.7

        # Start at the leftmost posture
        theta_tracker.set_value(theta_min)
        self.wait(0.5)

        # Sweep through the central textbook posture (θ2 = 0) to the rightmost posture
        self.play(theta_tracker.animate.set_value(0.0), run_time=3, rate_func=smooth)
        self.wait(0.5)
        self.play(theta_tracker.animate.set_value(theta_max), run_time=5, rate_func=linear)

        self.wait(2)


class FourBarAccelerationsFig410(MovingCameraScene):
    def construct(self):
        # Increase scene scale and stroke sizes for better visibility
        SCALE = 0.12
        VEL_SCALE = 0.02
        ACC_SCALE = 0.0016
        inch = 1.0

        R_BA = 4 * inch
        R_CB = 18 * inch
        R_CD = 11 * inch
        R_DA = 10 * inch
        R_GB = 10 * inch
        R_EC = 4 * inch
        R_HD = 7 * inch
        R_FH = 3 * inch

        omega2 = 94.25
        alpha2 = 0.0

        A = np.array([0.0, 0.0])
        D = np.array([R_DA, 0.0])
        theta2_0 = np.deg2rad(120)
        B = A + R_BA * np.array([np.cos(theta2_0), np.sin(theta2_0)])

        # Solve C by intersecting circle(D, R_CD) and circle(B, R_CB)
        def circle_intersections(p0, r0, p1, r1):
            d = np.linalg.norm(p1 - p0)
            if d < 1e-9:
                return []
            a = (r0**2 - r1**2 + d**2) / (2*d)
            h2 = r0**2 - a**2
            if h2 < 0:
                return []
            h = np.sqrt(max(h2, 0.0))
            p2 = p0 + a * (p1 - p0) / d
            rx = (p1[0] - p0[0]) / d
            ry = (p1[1] - p0[1]) / d
            offset = np.array([-ry * h, rx * h])
            return [p2 + offset, p2 - offset]

        candidates = circle_intersections(D, R_CD, B, R_CB)
        C = max(candidates, key=lambda p: p[1]) if candidates else D + np.array([R_CD, 0.0])

        # Unit directions along coupler and link 4
        dir_CB = (C - B) / (np.linalg.norm(C - B) + 1e-9)
        dir_DC = (C - D) / (np.linalg.norm(C - D) + 1e-9)

        # Points along specified distances
        E = C - R_EC * dir_CB
        G = B + R_GB * dir_CB
        H = D + R_HD * dir_DC
        # Make FH perpendicular to link 4, pointing downward
        perp4 = np.array([-dir_DC[1], dir_DC[0]])
        if perp4[1] > 0:
            perp4 = -perp4
        F = H + R_FH * perp4

        def to_point(p):
            return np.array([p[0], p[1], 0.0]) * SCALE

        def arrow(o, v, scale, color):
            return Arrow(to_point(o), to_point(o + v * scale), buff=0, stroke_width=3, color=color)

        ground = Line(to_point(np.array([-10.0, -2.0])), to_point(np.array([22.0, -2.0])), color=GREY_B, stroke_width=6)

        # Larger joint markers
        A_dot = Dot(to_point(A), color=WHITE, radius=0.08)
        B_dot = Dot(to_point(B), color=YELLOW, radius=0.09)
        C_dot = Dot(to_point(C), color=WHITE, radius=0.08)
        D_dot = Dot(to_point(D), color=WHITE, radius=0.08)
        E_dot = Dot(to_point(E), color=TEAL, radius=0.08)
        F_dot = Dot(to_point(F), color=ORANGE, radius=0.08)
        G_dot = Dot(to_point(G), color=BLUE, radius=0.08)
        H_dot = Dot(to_point(H), color=PURPLE, radius=0.08)

        # Thicker links
        link2 = Line(to_point(A), to_point(B), color=BLUE, stroke_width=8)
        coupler = Line(to_point(B), to_point(C), color=GREEN, stroke_width=8)
        link4 = Line(to_point(D), to_point(C), color=PURPLE, stroke_width=8)

        A_lbl = MathTex("A", font_size=28).next_to(A_dot, DOWN + LEFT, buff=0.1)
        B_lbl = MathTex("B", font_size=28).next_to(B_dot, UP + LEFT, buff=0.1)
        C_lbl = MathTex("C", font_size=28).next_to(C_dot, UP + RIGHT, buff=0.1)
        D_lbl = MathTex("D", font_size=28).next_to(D_dot, DOWN + RIGHT, buff=0.1)
        E_lbl = MathTex("E", font_size=26).next_to(E_dot, DOWN + RIGHT, buff=0.1)
        F_lbl = MathTex("F", font_size=26).next_to(F_dot, DOWN + RIGHT, buff=0.1)
        G_lbl = MathTex("G", font_size=26).next_to(G_dot, UP + LEFT, buff=0.1)
        H_lbl = MathTex("H", font_size=26).next_to(H_dot, DOWN + LEFT, buff=0.1)

        vA = omega2 * R_BA * np.array([-np.sin(theta2_0), np.cos(theta2_0)])
        aA = alpha2 * R_BA * np.array([-np.sin(theta2_0), np.cos(theta2_0)]) + (omega2**2) * R_BA * np.array([-np.cos(theta2_0), -np.sin(theta2_0)])

        vB = vA
        aB = aA

        omega3 = omega2 * 0.35
        alpha3 = 0.0
        omega4 = omega2 * 0.15
        alpha4 = 0.0

        theta3 = np.arctan2((C - B)[1], (C - B)[0])
        theta4 = np.arctan2((C - D)[1], (C - D)[0])

        vC = omega4 * np.linalg.norm(C - D) * np.array([-np.sin(theta4), np.cos(theta4)])
        aC = alpha4 * np.linalg.norm(C - D) * np.array([-np.sin(theta4), np.cos(theta4)]) + (omega4**2) * np.linalg.norm(C - D) * np.array([-np.cos(theta4), -np.sin(theta4)])

        vE = vC + omega3 * np.linalg.norm(C - E) * np.array([-np.sin(theta3), np.cos(theta3)])
        aE = aC + alpha3 * np.linalg.norm(C - E) * np.array([-np.sin(theta3), np.cos(theta3)]) + (omega3**2) * np.linalg.norm(C - E) * np.array([-np.cos(theta3), -np.sin(theta3)])

        vF = vC + omega3 * np.linalg.norm(C - F) * np.array([-np.sin(theta3), np.cos(theta3)])
        aF = aC + alpha3 * np.linalg.norm(C - F) * np.array([-np.sin(theta3), np.cos(theta3)]) + (omega3**2) * np.linalg.norm(C - F) * np.array([-np.cos(theta3), -np.sin(theta3)])

        velE = arrow(E, vE, VEL_SCALE, BLUE)
        accE = arrow(E, aE, ACC_SCALE, GREEN)
        velF = arrow(F, vF, VEL_SCALE, BLUE)
        accF = arrow(F, aF, ACC_SCALE, GREEN)

        omega_text = VGroup(
            MathTex("\\omega_2 = 94.25\\,\\text{rad/s}\\,\\text{ccw}", font_size=26),
            Text("alpha_3, alpha_4 shown qualitatively", font_size=24),
        ).arrange(DOWN).to_corner(UR)

        # Zoom in tighter on the mechanism
        self.camera.frame.set_width(6)
        # Center frame roughly around the coupler mid-point
        mid = (to_point(B) + to_point(C)) / 2
        self.camera.frame.move_to(mid)
        self.add(
            ground,
            link2, coupler, link4,
            A_dot, B_dot, C_dot, D_dot, E_dot, F_dot, G_dot, H_dot,
            A_lbl, B_lbl, C_lbl, D_lbl, E_lbl, F_lbl, G_lbl, H_lbl,
            velE, accE, velF, accF,
            omega_text,
        )
        self.wait(3)

class MechanismAnimation(Scene):
    def construct(self):
        # Scale factor (Manim units per inch)
        scale = 0.5
        
        # Given parameters (in inches, converted to Manim units)
        R_O2O = np.array([7.0, 5.0, 0]) * scale
        R_AO2 = 5.0 * scale
        R_BA = 2.75 * scale
        R_CD = 6.0 * scale
        R_BD = 2.818 * scale
        
        # Block dimensions (2 in x 1 in)
        block_width = 2.0 * scale
        block_height = 1.0 * scale
        
        # Fixed positions
        O2_pos = np.array([-4, -1, 0])
        O_pos = O2_pos + R_O2O
        
        # Ground hatching
        ground_line = Line(LEFT * 7 + DOWN * 3, RIGHT * 7 + DOWN * 3, color=GRAY)
        hatch_lines = VGroup(*[
            Line(LEFT * 7 + DOWN * 3, LEFT * 7 + DOWN * 3.3, color=GRAY),
            Line(LEFT * 6 + DOWN * 3, LEFT * 6 + DOWN * 3.3, color=GRAY),
            Line(LEFT * 5 + DOWN * 3, LEFT * 5 + DOWN * 3.3, color=GRAY),
            Line(LEFT * 4 + DOWN * 3, LEFT * 4 + DOWN * 3.3, color=GRAY),
            Line(LEFT * 3 + DOWN * 3, LEFT * 3 + DOWN * 3.3, color=GRAY),
            Line(LEFT * 2 + DOWN * 3, LEFT * 2 + DOWN * 3.3, color=GRAY),
            Line(LEFT * 1 + DOWN * 3, LEFT * 1 + DOWN * 3.3, color=GRAY),
            Line(ORIGIN + DOWN * 3, ORIGIN + DOWN * 3.3, color=GRAY),
            Line(RIGHT * 1 + DOWN * 3, RIGHT * 1 + DOWN * 3.3, color=GRAY),
            Line(RIGHT * 2 + DOWN * 3, RIGHT * 2 + DOWN * 3.3, color=GRAY),
            Line(RIGHT * 3 + DOWN * 3, RIGHT * 3 + DOWN * 3.3, color=GRAY),
            Line(RIGHT * 4 + DOWN * 3, RIGHT * 4 + DOWN * 3.3, color=GRAY),
            Line(RIGHT * 5 + DOWN * 3, RIGHT * 5 + DOWN * 3.3, color=GRAY),
            Line(RIGHT * 6 + DOWN * 3, RIGHT * 6 + DOWN * 3.3, color=GRAY),
            Line(RIGHT * 7 + DOWN * 3, RIGHT * 7 + DOWN * 3.3, color=GRAY),
        ])
        
        # Fixed pivots
        O2_dot = Dot(O2_pos, color=WHITE, radius=0.1)
        O2_label = Text("O₂", font_size=24).next_to(O2_dot, DOWN)
        O2_ground = DashedLine(O2_pos, O2_pos + DOWN * 2, color=GRAY)
        
        O_dot = Dot(O_pos, color=WHITE, radius=0.1)
        O_label = Text("O", font_size=24).next_to(O_dot, DOWN)
        
        # Slider rail (horizontal)
        slider_y = 0.5
        slider_rail = Line(LEFT * 6 + UP * slider_y, RIGHT * 6 + UP * slider_y, 
                          color=GRAY, stroke_width=6)
        
        # Create labels for dimensions
        title = Text("Slider-Crank Mechanism", font_size=36).to_edge(UP)
        params = VGroup(
            Text("R_AO₂ = 5.0 in", font_size=20),
            Text("R_BA = 2.75 in", font_size=20),
            Text("R_BD = 2.818 in", font_size=20),
            Text("R_CD = 6.0 in", font_size=20),
        ).arrange(DOWN, aligned_edge=LEFT).to_corner(UL, buff=0.5).shift(DOWN * 0.5)
        
        self.add(ground_line, hatch_lines)
        self.add(O2_dot, O2_label, O2_ground, O_dot, O_label)
        self.add(slider_rail)
        self.add(title, params)
        
        # Animation function
        def get_mechanism_positions(theta):
            # Point A rotates around O2
            A_pos = O2_pos + R_AO2 * np.array([np.cos(theta), np.sin(theta), 0])
            
            # Solve for B position using geometric constraints
            # B is R_BA from A and R_BD from slider rail
            # Using circle-circle intersection
            
            # B must be on a circle of radius R_BA centered at A
            # B must also be R_BD from the slider rail (vertical distance)
            
            # Solve for B_x given constraints
            # Distance from A to B = R_BA
            # |B_y - slider_y| ≈ R_BD (approximately, for simplified kinematics)
            
            # Using iterative approach for B position
            B_y = slider_y + R_BD * 0.8  # Approximate vertical position
            dx_squared = R_BA**2 - (B_y - A_pos[1])**2
            
            if dx_squared > 0:
                B_x = A_pos[0] + np.sqrt(dx_squared)
            else:
                B_x = A_pos[0] + R_BA
            
            B_pos = np.array([B_x, B_y, 0])
            
            # D is on the slider rail, vertically aligned with B
            D_pos = np.array([B_x, slider_y, 0])
            
            # C is R_CD to the right of D on the slider
            C_pos = D_pos + np.array([R_CD, 0, 0])
            
            return A_pos, B_pos, D_pos, C_pos
        
        # Initial positions
        theta_tracker = ValueTracker(PI/4)
        
        A_pos, B_pos, D_pos, C_pos = get_mechanism_positions(theta_tracker.get_value())
        
        # Create mechanism components
        # Link O2-A (crank)
        link_O2A = always_redraw(lambda: Line(
            O2_pos, 
            get_mechanism_positions(theta_tracker.get_value())[0],
            color=RED, stroke_width=8
        ))
        
        A_dot = always_redraw(lambda: Dot(
            get_mechanism_positions(theta_tracker.get_value())[0],
            color=RED, radius=0.08
        ))
        
        A_label = always_redraw(lambda: Text("A", font_size=20).next_to(
            get_mechanism_positions(theta_tracker.get_value())[0], UP
        ))
        
        # Link A-B (coupler)
        link_AB = always_redraw(lambda: Line(
            get_mechanism_positions(theta_tracker.get_value())[0],
            get_mechanism_positions(theta_tracker.get_value())[1],
            color=BLUE, stroke_width=8
        ))
        
        B_dot = always_redraw(lambda: Dot(
            get_mechanism_positions(theta_tracker.get_value())[1],
            color=BLUE, radius=0.08
        ))
        
        B_label = always_redraw(lambda: Text("B", font_size=20).next_to(
            get_mechanism_positions(theta_tracker.get_value())[1], UP
        ))
        
        # Block 5 (at B)
        block_5 = always_redraw(lambda: Rectangle(
            width=block_width, height=block_height,
            color=GRAY, fill_opacity=0.5, stroke_width=2
        ).move_to(get_mechanism_positions(theta_tracker.get_value())[1]))
        
        # Link B-D
        link_BD = always_redraw(lambda: Line(
            get_mechanism_positions(theta_tracker.get_value())[1],
            get_mechanism_positions(theta_tracker.get_value())[2],
            color=PURPLE, stroke_width=8
        ))
        
        D_dot = always_redraw(lambda: Dot(
            get_mechanism_positions(theta_tracker.get_value())[2],
            color=PURPLE, radius=0.08
        ))
        
        D_label = always_redraw(lambda: Text("D", font_size=20).next_to(
            get_mechanism_positions(theta_tracker.get_value())[2], DOWN
        ))
        
        # Block 6 (at D)
        block_6 = always_redraw(lambda: Rectangle(
            width=block_width, height=block_height,
            color=GRAY, fill_opacity=0.5, stroke_width=2
        ).move_to(get_mechanism_positions(theta_tracker.get_value())[2]))
        
        # Link D-C
        link_DC = always_redraw(lambda: Line(
            get_mechanism_positions(theta_tracker.get_value())[2],
            get_mechanism_positions(theta_tracker.get_value())[3],
            color=GREEN, stroke_width=8
        ))
        
        C_dot = always_redraw(lambda: Dot(
            get_mechanism_positions(theta_tracker.get_value())[3],
            color=GREEN, radius=0.08
        ))
        
        C_label = always_redraw(lambda: Text("C", font_size=20).next_to(
            get_mechanism_positions(theta_tracker.get_value())[3], UP
        ))
        
        # Add all components
        self.add(link_O2A, link_AB, link_BD, link_DC)
        self.add(block_5, block_6)
        self.add(A_dot, B_dot, D_dot, C_dot)
        self.add(A_label, B_label, D_label, C_label)
        
        # Animate the mechanism
        self.play(
            theta_tracker.animate.set_value(theta_tracker.get_value() + 2*PI),
            rate_func=linear,
            run_time=8
        )
        
        self.wait(1)


from manim import *
import numpy as np

class FourBarLinkage(Scene):
    def construct(self):
        # Given parameters
        omega_2 = 900  # rev/min ccw
        omega_2_rad = omega_2 * 2 * PI / 60  # Convert to rad/s (94.25 rad/s)
        
        # Scale factor for visualization (reduced)
        scale = 1.2
        
        # Link lengths (estimated from figure)
        L_AB = 2.0 * scale  # Crank length (link 2)
        L_BC = 5.0 * scale  # Coupler length (link 3)
        L_CD = 3.5 * scale  # Output link length (link 4)
        L_AD = 6.0 * scale  # Ground link length
        
        # Perpendicular distances
        L_GE = 2.0 * scale  # E is perpendicular to BC at G
        L_FH = 1.5 * scale  # F is perpendicular to CD at H
        
        # Initial angle of crank (120° from horizontal)
        theta_2_initial = 120 * DEGREES
        
        # Time tracker
        time_tracker = ValueTracker(0)
        
        # Fixed ground points - A at origin BEFORE axes shift
        A_pos = ORIGIN
        D_pos = A_pos + np.array([L_AD, 0, 0])
        
        # Create coordinate system centered at A
        axes = Axes(
            x_range=[-1, 9, 1],
            y_range=[-3, 4, 1],
            x_length=8,
            y_length=5,
            axis_config={"color": GRAY, "include_tip": True},
        )
        
        # Shift everything down to prevent overflow
        shift_amount = DOWN * 1.2 + LEFT * 0.5
        axes.shift(shift_amount)
        A_pos = A_pos + shift_amount
        D_pos = D_pos + shift_amount
        
        x_label = MathTex("x_1").next_to(axes.x_axis.get_end(), RIGHT)
        y_label = MathTex("y_1").next_to(axes.y_axis.get_end(), UP)
        
        # Ground hatching for fixed pivots
        def create_ground_hatch(pos):
            hatch = VGroup()
            for i in range(-5, 6):
                line = Line(
                    pos + np.array([i * 0.15, -0.3, 0]),
                    pos + np.array([i * 0.15 + 0.2, -0.5, 0]),
                    color=GRAY, stroke_width=2
                )
                hatch.add(line)
            ground_line = Line(pos + LEFT * 0.8 + DOWN * 0.3, 
                             pos + RIGHT * 0.8 + DOWN * 0.3, 
                            color=GRAY, stroke_width=4)
            return VGroup(hatch, ground_line)
        
        ground_A = create_ground_hatch(A_pos)
        ground_D = create_ground_hatch(D_pos)
        
        # Fixed pivot points
        A_dot = Dot(A_pos, color=WHITE, radius=0.12)
        A_label = MathTex("A", font_size=36).next_to(A_dot, DOWN + LEFT, buff=0.1)
        
        D_dot = Dot(D_pos, color=WHITE, radius=0.12)
        D_label = MathTex("D", font_size=36).next_to(D_dot, DOWN)
        
        # Title and parameters
        title = Text("Four-Bar Linkage Mechanism", font_size=40).to_edge(UP)
        params = VGroup(
            MathTex(r"\omega_2 = 900 \text{ rev/min ccw}", font_size=28),
            MathTex(r"\omega_2 = 94.25 \text{ rad/s}", font_size=28),
        ).arrange(DOWN, aligned_edge=LEFT).to_corner(UR, buff=0.5)
        
        # Live angle display
        angle_display = always_redraw(lambda: MathTex(
            r"\theta_2 = " + f"{(theta_2_initial + omega_2_rad * time_tracker.get_value()) * 180 / PI % 360:.1f}" + r"°",
            font_size=32,
            color=YELLOW
        ).next_to(params, DOWN, aligned_edge=LEFT, buff=0.3))
        
        # Time display
        time_display = always_redraw(lambda: MathTex(
            r"t = " + f"{time_tracker.get_value():.2f}" + r" \text{ s}",
            font_size=28,
            color=GREEN
        ).next_to(angle_display, DOWN, aligned_edge=LEFT, buff=0.2))
        
        # Function to calculate linkage positions
        def get_linkage_positions(t):
            theta_2 = theta_2_initial + omega_2_rad * t
            
            # Point B (end of crank)
            B_pos = A_pos + L_AB * np.array([np.cos(theta_2), np.sin(theta_2), 0])
            
            # Point C using geometric constraint (intersection of two circles)
            dx = D_pos[0] - B_pos[0]
            dy = D_pos[1] - B_pos[1]
            d = np.sqrt(dx**2 + dy**2)
            
            # Check if solution exists
            if d > L_BC + L_CD or d < abs(L_BC - L_CD):
                # Use approximate solution
                angle_BD = np.arctan2(dy, dx)
                C_pos = B_pos + L_BC * np.array([np.cos(angle_BD), np.sin(angle_BD), 0])
            else:
                # Solve for C position
                a = (L_BC**2 - L_CD**2 + d**2) / (2 * d)
                h = np.sqrt(max(0, L_BC**2 - a**2))
                
                P_x = B_pos[0] + a * (dx) / d
                P_y = B_pos[1] + a * (dy) / d
                
                C_pos = np.array([
                    P_x + h * (-dy) / d,
                    P_y + h * (dx) / d,
                    0
                ])
            
            # Point G (on link BC, about 60% from B to C)
            G_pos = B_pos + 0.6 * (C_pos - B_pos)
            
            # Point E: perpendicular to BC at G (DOWNWARD - negative perpendicular)
            # Direction of BC
            BC_dir = (C_pos - B_pos) / np.linalg.norm(C_pos - B_pos)
            # Perpendicular direction (rotate 90° CW for downward)
            perp_dir = np.array([BC_dir[1], -BC_dir[0], 0])
            E_pos = G_pos + L_GE * perp_dir
            
            # Point H: midpoint of CD
            H_pos = (C_pos + D_pos) / 2
            
            # Point F: perpendicular to CD at H
            # Direction of CD
            CD_dir = (D_pos - C_pos) / np.linalg.norm(D_pos - C_pos)
            # Perpendicular direction (rotate 90° CCW)
            perp_dir_CD = np.array([-CD_dir[1], CD_dir[0], 0])
            F_pos = H_pos + L_FH * perp_dir_CD
            
            return B_pos, C_pos, E_pos, G_pos, F_pos, H_pos, theta_2
        
        # Create linkage elements with always_redraw
        # Link 2 (Crank AB) - Red
        link_2 = always_redraw(lambda: Line(
            A_pos,
            get_linkage_positions(time_tracker.get_value())[0],
            color=RED, stroke_width=10
        ))
        
        # Link 3 (Coupler BC) - Blue with fill
        def create_coupler():
            B, C, E, G, F, H, theta = get_linkage_positions(time_tracker.get_value())
            coupler = Polygon(
                B, C,
                color=BLUE, fill_opacity=0.3, stroke_width=8
            )
            return coupler
        
        link_3 = always_redraw(create_coupler)
        
        # Link 4 (Output CD) - Blue 
        link_4 = always_redraw(lambda: Line(
            get_linkage_positions(time_tracker.get_value())[1],
            D_pos,
            color=BLUE_C, stroke_width=10
        ))
        
        # GE perpendicular line (dashed)
        line_GE = always_redraw(lambda: DashedLine(
            get_linkage_positions(time_tracker.get_value())[3],
            get_linkage_positions(time_tracker.get_value())[2],
            color=YELLOW, stroke_width=3
        ))
        
        # FH perpendicular line (dashed)
        line_FH = always_redraw(lambda: DashedLine(
            get_linkage_positions(time_tracker.get_value())[5],
            get_linkage_positions(time_tracker.get_value())[4],
            color=GREEN, stroke_width=3
        ))
        
        # Joint points
        B_dot = always_redraw(lambda: Dot(
            get_linkage_positions(time_tracker.get_value())[0],
            color=WHITE, radius=0.1
        ))
        
        B_label = always_redraw(lambda: MathTex("B", font_size=32).next_to(
            get_linkage_positions(time_tracker.get_value())[0], LEFT
        ))
        
        C_dot = always_redraw(lambda: Dot(
            get_linkage_positions(time_tracker.get_value())[1],
            color=WHITE, radius=0.1
        ))
        
        C_label = always_redraw(lambda: MathTex("C", font_size=32).next_to(
            get_linkage_positions(time_tracker.get_value())[1], UP + RIGHT
        ))
        
        # Point E
        E_dot = always_redraw(lambda: Dot(
            get_linkage_positions(time_tracker.get_value())[2],
            color=YELLOW, radius=0.08
        ))
        
        E_label = always_redraw(lambda: MathTex("E", font_size=28, color=YELLOW).next_to(
            get_linkage_positions(time_tracker.get_value())[2], DOWN
        ))
        
        # Point F
        F_dot = always_redraw(lambda: Dot(
            get_linkage_positions(time_tracker.get_value())[4],
            color=GREEN, radius=0.08
        ))
        
        F_label = always_redraw(lambda: MathTex("F", font_size=28, color=GREEN).next_to(
            get_linkage_positions(time_tracker.get_value())[4], RIGHT
        ))
        
        # Point G
        G_dot = always_redraw(lambda: Dot(
            get_linkage_positions(time_tracker.get_value())[3],
            color=ORANGE, radius=0.08
        ))
        
        G_label = always_redraw(lambda: MathTex("G", font_size=28, color=ORANGE).next_to(
            get_linkage_positions(time_tracker.get_value())[3], DOWN + LEFT, buff=0.1
        ))
        
        # Point H
        H_dot = always_redraw(lambda: Dot(
            get_linkage_positions(time_tracker.get_value())[5],
            color=PURPLE, radius=0.08
        ))
        
        H_label = always_redraw(lambda: MathTex("H", font_size=28, color=PURPLE).next_to(
            get_linkage_positions(time_tracker.get_value())[5], LEFT
        ))
        
        # Angle annotation at A with arrow (dynamic)
        def create_angle_arc_with_arrow():
            current_theta = get_linkage_positions(time_tracker.get_value())[6]
            arc_radius = 0.7
            
            # Create arc with tip (built-in arrow functionality)
            arc = ArcBetweenPoints(
                A_pos + arc_radius * RIGHT,
                A_pos + arc_radius * np.array([np.cos(current_theta), np.sin(current_theta), 0]),
                angle=current_theta,
                color=YELLOW,
                stroke_width=3
            )
            
            # Calculate tangent direction at the end of arc
            tangent_angle = current_theta + PI/2  # Perpendicular to radius
            arrow_direction = np.array([np.cos(tangent_angle), np.sin(tangent_angle), 0])
            
            # Create arrow tip at end of arc
            end_point = A_pos + arc_radius * np.array([np.cos(current_theta), np.sin(current_theta), 0])
            
            arrow = Arrow(
                end_point - 0.001 * arrow_direction,
                end_point + 0.25 * arrow_direction,
                color=YELLOW,
                buff=0,
                stroke_width=3,
                tip_length=0.2,
                max_stroke_width_to_length_ratio=10,
                max_tip_length_to_length_ratio=1
            )
            
            return VGroup(arc, arrow)
        
        angle_arc = always_redraw(create_angle_arc_with_arrow)
        
        # Omega label near the angle arc
        omega_label = always_redraw(lambda: MathTex(
            r"\omega_2", font_size=24, color=YELLOW
        ).move_to(A_pos + np.array([1.0, 0.6, 0])))
        
        # Link number labels
        link_2_num = always_redraw(lambda: MathTex("2", font_size=28, color=RED).move_to(
            (A_pos + get_linkage_positions(time_tracker.get_value())[0]) / 2 + DOWN * 0.4
        ))
        
        link_3_num = always_redraw(lambda: MathTex("3", font_size=28, color=BLUE).move_to(
            (get_linkage_positions(time_tracker.get_value())[0] + 
             get_linkage_positions(time_tracker.get_value())[1]) / 2 + UP * 0.3
        ))
        
        link_4_num = always_redraw(lambda: MathTex("4", font_size=28, color=BLUE_C).move_to(
            (D_pos + get_linkage_positions(time_tracker.get_value())[1]) / 2 + RIGHT * 0.5
        ))
        
        # Right angle markers for perpendicularity
        def create_right_angle_marker(corner, p1, p2, size=0.2):
            v1 = (p1 - corner) / np.linalg.norm(p1 - corner) * size
            v2 = (p2 - corner) / np.linalg.norm(p2 - corner) * size
            marker = VGroup(
                Line(corner + v1, corner + v1 + v2, stroke_width=2),
                Line(corner + v1 + v2, corner + v2, stroke_width=2)
            )
            return marker
        
        right_angle_G = always_redraw(lambda: create_right_angle_marker(
            get_linkage_positions(time_tracker.get_value())[3],  # G
            get_linkage_positions(time_tracker.get_value())[0],  # B
            get_linkage_positions(time_tracker.get_value())[2],  # E
        ).set_color(YELLOW))
        
        right_angle_H = always_redraw(lambda: create_right_angle_marker(
            get_linkage_positions(time_tracker.get_value())[5],  # H
            get_linkage_positions(time_tracker.get_value())[1],  # C
            get_linkage_positions(time_tracker.get_value())[4],  # F
        ).set_color(GREEN))
        
        # Add all elements to scene
        self.add(axes, x_label, y_label)
        self.add(ground_A, ground_D)
        self.add(A_dot, A_label, D_dot, D_label)
        self.add(title, params, angle_display, time_display)
        
        # Add links (order matters for layering)
        self.add(link_4, link_3, link_2)
        self.add(line_GE, line_FH)
        
        # Add points and labels
        self.add(B_dot, C_dot, E_dot, F_dot, G_dot, H_dot)
        self.add(B_label, C_label, E_label, F_label, G_label, H_label)
        self.add(angle_arc, omega_label)
        self.add(link_2_num, link_3_num, link_4_num)
        self.add(right_angle_G, right_angle_H)
        
        # Animate the mechanism - two full rotations
        rotation_time = (4 * PI) / omega_2_rad  # Time for 2 full rotations
        
        self.play(
            time_tracker.animate.set_value(rotation_time),
            rate_func=linear,
            run_time=10
        )
        
        self.wait(1)