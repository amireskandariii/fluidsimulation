import numpy as np
from vispy import app, gloo

# Simulation parameters
N = 200  # Resolution of the grid (NxN cells)
dt = 10  # Time step for the simulation
viscosity = 0.001  # Viscosity of the fluid
diffusion = 0.0  # Diffusion rate for the density field

# Initialize fields
# Velocity fields (u, v) and their previous states (u_prev, v_prev)
u = np.zeros((N+2, N+2), dtype=np.float32)
v = np.zeros((N+2, N+2), dtype=np.float32)
u_prev = np.zeros((N+2, N+2), dtype=np.float32)
v_prev = np.zeros((N+2, N+2), dtype=np.float32)
# Density field and its previous state
density = np.zeros((N+2, N+2), dtype=np.float32)
density_prev = np.zeros((N+2, N+2), dtype=np.float32)

# Fluid solver functions
# Function to add source field to another field
def add_source(x, s):
    x += dt * s

# Function to set boundary conditions for velocity and density fields
def set_bnd(b, x):
    for i in range(1, N+1):
        # Set boundary values for horizontal edges (b==1)
        x[0, i] = -x[1, i] if b == 1 else x[1, i]
        x[N+1, i] = -x[N, i] if b == 1 else x[N, i]
        # Set boundary values for vertical edges (b==2)
        x[i, 0] = -x[i, 1] if b == 2 else x[i, 1]
        x[i, N+1] = -x[i, N] if b == 2 else x[i, N]
    # Set corner values
    x[0, 0] = 0.5 * (x[1, 0] + x[0, 1])
    x[0, N+1] = 0.5 * (x[1, N+1] + x[0, N])
    x[N+1, 0] = 0.5 * (x[N, 0] + x[N+1, 1])
    x[N+1, N+1] = 0.5 * (x[N, N+1] + x[N+1, N])

# Diffusion step: Applies diffusion to a field
def diffuse(b, x, x0, diff):
    a = dt * diff * N * N  # Diffusion coefficient
    for _ in range(20):  # Iterative solver for diffusion
        x[1:-1, 1:-1] = (x0[1:-1, 1:-1] + a * (x[2:, 1:-1] + x[:-2, 1:-1] +
                                                x[1:-1, 2:] + x[1:-1, :-2])) / (1 + 4 * a)
        set_bnd(b, x)  # Set boundary conditions

# Advection step: Moves the fluid based on velocity fields
def advect(b, d, d0, u, v):
    dt0 = dt * N  # Scaling factor for velocity
    for i in range(1, N+1):
        for j in range(1, N+1):
            # Trace back to find the source cell contributing to current cell
            x = i - dt0 * u[i, j]
            y = j - dt0 * v[i, j]
            # Clamp coordinates to avoid going out of bounds
            if x < 0.5:
                x = 0.5
            if x > N + 0.5:
                x = N + 0.5
            i0 = int(x)
            i1 = i0 + 1
            if y < 0.5:
                y = 0.5
            if y > N + 0.5:
                y = N + 0.5
            j0 = int(y)
            j1 = j0 + 1
            # Linear interpolation to determine the value at the current cell
            s1 = x - i0
            s0 = 1 - s1
            t1 = y - j0
            t0 = 1 - t1
            d[i, j] = (s0 * (t0 * d0[i0, j0] + t1 * d0[i0, j1]) +
                       s1 * (t0 * d0[i1, j0] + t1 * d0[i1, j1]))
    set_bnd(b, d)  # Set boundary conditions

# Projection step: Enforces mass conservation (incompressibility)
def project(u, v, p, div):
    h = 1.0 / N
    # Compute the divergence of the velocity field
    div[1:-1, 1:-1] = -0.5 * h * (u[2:, 1:-1] - u[:-2, 1:-1] +
                                  v[1:-1, 2:] - v[1:-1, :-2])
    p[1:-1, 1:-1] = 0  # Initialize the pressure field
    set_bnd(0, div)
    set_bnd(0, p)
    # Solve for pressure
    for _ in range(20):
        p[1:-1, 1:-1] = (div[1:-1, 1:-1] + p[2:, 1:-1] + p[:-2, 1:-1] +
                         p[1:-1, 2:] + p[1:-1, :-2]) / 4
        set_bnd(0, p)
    # Subtract pressure gradient from velocity to make it divergence-free
    u[1:-1, 1:-1] -= 0.5 * (p[2:, 1:-1] - p[:-2, 1:-1]) / h
    v[1:-1, 1:-1] -= 0.5 * (p[1:-1, 2:] - p[1:-1, :-2]) / h
    set_bnd(1, u)
    set_bnd(2, v)

# Function to perform a velocity step (diffusion and advection for velocity)
def vel_step(u, v, u0, v0):
    add_source(u, u0)  # Add velocity sources
    add_source(v, v0)

    u0[:], v0[:] = u.copy(), v.copy()
    diffuse(1, u, u0, viscosity)  # Diffuse velocity
    diffuse(2, v, v0, viscosity)

    project(u, v, u0, v0)  # Enforce incompressibility

    u0[:], v0[:] = u.copy(), v.copy()
    advect(1, u, u0, u0, v0)  # Advect velocity
    advect(2, v, v0, u0, v0)

    project(u, v, u0, v0)  # Enforce incompressibility again

# Function to perform a density step (diffusion and advection for density)
def dens_step(x, x0, u, v):
    add_source(x, x0)  # Add density sources

    x0[:] = x.copy()
    diffuse(0, x, x0, diffusion)  # Diffuse density

    x0[:] = x.copy()
    advect(0, x, x0, u, v)  # Advect density

# Function to add density to a specific position
def add_density(x_pos, y_pos, amount):
    density_prev[x_pos, y_pos] += amount

# Function to add velocity to a specific position
def add_velocity(x_pos, y_pos, amount_u, amount_v):
    u_prev[x_pos, y_pos] += amount_u
    v_prev[x_pos, y_pos] += amount_v

# Visualization code using VisPy
# Vertex shader for rendering the density field
vertex_shader = """
attribute vec2 a_position;
varying vec2 v_texcoord;
void main() {
    v_texcoord = a_position;
    gl_Position = vec4(2.0 * a_position - 1.0, 0.0, 1.0);
}
"""

# Fragment shader for rendering the density field
fragment_shader = """
uniform sampler2D u_texture;
varying vec2 v_texcoord;
void main() {
    float d = texture2D(u_texture, v_texcoord).r;
    gl_FragColor = vec4(d, d, d, 1.0);  // Gray scale based on density
}
"""

# Canvas class to handle visualization
class Canvas(app.Canvas):
    def __init__(self):
        app.Canvas.__init__(self, keys='interactive', size=(600, 600))
        # Compile shaders and set up GPU program
        self.program = gloo.Program(vertex_shader, fragment_shader)
        self.texture = gloo.Texture2D((N, N), format='luminance', interpolation='linear')
        self.program['u_texture'] = self.texture

        # Full-screen quad covering the viewport
        positions = np.array([
            [0.0, 0.0],  # Bottom-left
            [1.0, 0.0],  # Bottom-right
            [0.0, 1.0],  # Top-left
            [1.0, 1.0],  # Top-right
        ], dtype=np.float32)
        self.program['a_position'] = positions

        # Indices for drawing two triangles to form a quad
        self.indices = gloo.IndexBuffer([0, 1, 2, 1, 3, 2])

        # Set viewport and state for rendering
        gloo.set_viewport(0, 0, *self.physical_size)
        gloo.set_state(clear_color='white', blend=True, blend_func=('src_alpha', 'one_minus_src_alpha'))

        # Timer to trigger simulation updates
        self.timer = app.Timer('auto', connect=self.on_timer, start=True)

    def on_draw(self, event):
        gloo.clear()  # Clear the screen
        # Update texture data with current density field
        texture_data = np.flipud(density[1:-1, 1:-1])  # Flip data vertically
        self.texture.set_data(texture_data)
        self.program.draw('triangles', self.indices)  # Draw the quad

    def on_timer(self, event):
        # Simulation step
        global u, v, u_prev, v_prev, density, density_prev

        # Add sources (you can customize the positions and amounts)
        add_density(N//2, N//2, 100.0)  # Add density at the center
        add_velocity(N//2, N//2, np.random.uniform(-1, 1), np.random.uniform(-1, 1))  # Add random velocity at the center

        # Perform simulation steps for velocity and density
        vel_step(u, v, u_prev, v_prev)
        dens_step(density, density_prev, u, v)

        # Reset previous fields
        u_prev.fill(0)
        v_prev.fill(0)
        density_prev.fill(0)

        self.update()  # Request a redraw

# Run the application
if __name__ == '__main__':
    c = Canvas()
    c.show()
    app.run()
