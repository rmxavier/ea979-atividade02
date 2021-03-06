# Renders a 2D model into a PPM image
import sys
import numpy as np

# ---------- Configuration types and constants ----------

MAX_SIZE = 1024
MAX_VAL = 255
MAX_LINE_LEN = 10240-1 # 10240 characters minus the \0 terminator
DEFAULT_BACKGROUND = 255
CHANNELS_N = 3
COORD_N = 3
DEFAULT_COLOR = (0, 0, 0,)
IMAGE_DTYPE = np.uint8
VIEWPORT_DTYPE = np.int64
MODEL_DTYPE = np.float64
ZBUFFER_DTYPE = np.float64
ZBUFFER_BACKGROUND = -np.inf

# Geometries
GEOMETRY_REGION = "REGION"
GEOMETRY_POLYLINE = "POLYLINE"

# Commands
SET_BACKGROUND_COLOR = "c"
SET_COLOR = "C"
DRAW_GEOMETRY_LINE = "L"
DRAW_GEOMETRY_POLYLINE = "P"
DRAW_GEOMETRY_REGION = "R"
DRAW_GEOMETRY_SPHERE = "SPH"
DRAW_GEOMETRY_CUBE = "CUB"
SET_TRANSFORMATION_MATRIX = "M"
SET_VIEWPORT_MATRIX = "V"
APPLY_TRANSFORMATION_MATRIX = "m"
PUSH_TRANSFORMATION_MATRIX = "PUSH"
POP_TRANSFORMATION_MATRIX = "POP"


# ---------- Output routines ----------

def put_string(output, output_file):
    output = output.encode('ascii') if isinstance(output, str) else output
    written_n = output_file.write(output)
    if written_n != len(output):
        print('error writing to output stream', file=sys.stderr)
        sys.exit(1)


def save_ppm(image, output_file):
    # Defines image header
    magic_number_1 = 'P'
    magic_number_2 = '6'
    width  = image.shape[1]
    height = image.shape[0]
    end_of_header = '\n'

    # Writes header
    put_string(magic_number_1, output_file)
    put_string(magic_number_2, output_file)
    put_string('\n', output_file)
    put_string('%d %d\n' % (width, height), output_file)
    put_string('%d' % MAX_VAL, output_file)
    put_string(end_of_header, output_file)

    # Outputs image
    put_string(image.tobytes(), output_file)

# ---------- Auxiliar functions ----------
    
##
# Transform the input from a list of strings to a matrix of vectors in homogeneous coordinates
##
def get_points(p_list):
    p=[]
    s = [float(s) for s in p_list]
    # separates the input form of L comand, and R and P comand
    if len(p_list) % 3==0 :
        for i in range (0,2):
            p.append([s[3*i],s[3*i+1],s[3*i+2],1.0])
    else:
        for i in range (0,int(p_list[0])):
            p.append([s[3*i+1],s[3*i+2],s[3*i+3],1.0])
    return np.array(p)

##
# Transform the input from a list of strings to a matrix
##
def get_matrix(m_list):   
    m=[]
    s = [float(s) for s in m_list]
    for i in range (0,4):
        m.append([s[4*i],s[4*i+1],s[4*i+2],s[4*i+3]])
    return np.array(m)


##
# Applies a transformation matrix to a list of vectors
##   
def transform(matrix,vector):
    p = matrix.dot(vector.transpose())
    return p.transpose()

##
# Calculates the sphere points
##    
def get_sph(R,N_m,N_p):
    
    delta_teta = 180/(N_m)   
    delta_phi = 180/(N_p+1) 
    phi = delta_phi
    # Calculates the first circle , in xy plan
    l_c=([[0,R,0,1],[0,-R,0,1]])   
    for i in range(0,N_p):
        x = R*np.sin(np.radians(phi))
        y = R*np.cos(np.radians(phi))
        l_c.insert(-(i+1),[x,y,0,1])
        l_c.append([-x,-y,0,1])
        phi+=delta_phi
    # Apply a rotation matrix in axle y to get the other circle making the sphere   
    l_sph=np.array([l_c]) 
    a = np.cos(np.radians(delta_teta))
    b = np.sin(np.radians(delta_teta))
    M_r= np.array([[a,0,b,0],[0,1,0,0],[-b,0,a,0],[0,0,0,1]])
    for i in range (1,N_m):
        p = M_r.dot(l_sph[i-1].transpose())
        l_sph=np.append(l_sph,[p.transpose()], axis = 0) 
    return l_sph            # each position off this final array, is a matrix that carry the points of one circle

##
# Applies a transformation matrix for every circle that compound the sphere
##   
def transform_sph(matrix,array_sph):
    p = matrix.dot(array_sph[0].transpose())
    lt_sph = np.array([p.transpose()])
    for i in range(1,len(array_sph)):
        p = matrix.dot(array_sph[i].transpose())
        lt_sph = np.append(lt_sph,[p.transpose()], axis = 0) 
    return lt_sph

# ---------- Drawing/model routines ----------
def draw_line(image, V0, V1, color):
    #Transform the vectors in homogeneous coordinates to R3 int
    x0 = round((V0[0]/V0[3])+width/2)
    x1 = round((V1[0]/V1[3])+width/2)
    y0 = round((V0[1]/V0[3])+height/2)
    y1 = round((V1[1]/V1[3])+height/2)
    # Stores the z information
    z0 = round(V0[2])
    z1 = round(V1[2])
    dx = x1-x0
    dy = y1-y0
    dc = abs(dx) # delta x in book - here we are using row, col coordinates
    dr = abs(dy) # delta y in book
    if dr==0 and dc==0: #(evita casos err??nios que tivemos com os polos da esfera que ao arrendondar as coordenadas ficam as mesmas)
        return
    elif dr <= dc:
        # Line inclination is at most 1
        # Swaps points if c1<c0 and converts x,y coordinates to row,col coordinates
        # dx>=0 => x1>=x0 => c1>=x0
        r0 = height-1-y0 if dx>=0 else height-1-y1
        r1 = height-1-y1 if dx>=0 else height-1-y0
        c0 =     x0 if dx>=0 else x1
        c1 =     x1 if dx>=0 else x0
        # Implements Bresenham's midpoint algorithm for lines
        # (Klawonn. Introduction to Computer Graphics. 2nd Edition. Section 4.2, pp. 45???53)
        # ...deltas of Bressenham's algorithm
        d_horizontal = 2*dr      # delta east in book
        d_diagonal   = 2*(dr-dc) # delta northeast in book
        # ...draws line
        pixel_r = r0
        step_row = 1 if r1>=r0 else -1
        d = 2*dr - dc # starting D value, D_init in book
        for pixel_c in range(c0, c1+1):
            depth = z0 + (z1-z0)*(pixel_c-c0)/(c1-c0)
            # Check z-buffer conditions and if the pixel is inside the view port (solu????o para o clip longe do ideal, mas funciona)
            if pixel_c >=0 and pixel_c < width and pixel_r >=0 and pixel_r < height and depth>= zbuffer[pixel_r, pixel_c] :
                zbuffer[pixel_r, pixel_c] = depth
                image[pixel_r, pixel_c] = color
            if d<=0:
                d += d_horizontal
            else:
                d += d_diagonal
                pixel_r += step_row
    else:
        # Line inclination is greater than one -- inverts the roles of row and column
        # Swaps points if y1>y0 and converts x,y coordinates to row,col coordinates
        # dy<=0 => y1<=y0 => r1>=r0
        r0 = height-1-y0 if dy<=0 else height-1-y1
        r1 = height-1-y1 if dy<=0 else height-1-y0
        c0 =     x0 if dy<=0 else x1
        c1 =     x1 if dy<=0 else x0
        # Implements Bresenham's midpoint algorithm for lines
        # (Klawonn. Introduction to Computer Graphics. 2nd Edition. Section 4.2, pp. 45???53)
        # ...deltas of Bressenham's algorithm - same as above, but with coordinates inverted
        d_vertical = 2*dc
        d_diagonal = 2*(dc-dr)
        pixel_r = r0
        pixel_c = c0
        step_col = 1 if c1>=c0 else -1
        d = 2*dc - dr # starting D value, D_init in book
        for pixel_r in range(r0, r1+1):
            # Calculates the depth
            depth = z0 + (z1-z0)*(pixel_r-r0)/(r1-r0) 
            # Check z-buffer conditions and if the pixel is inside the view port   
            if pixel_c >=0 and pixel_c < width and pixel_r >=0 and pixel_r < height and depth>= zbuffer[pixel_r, pixel_c] :
                zbuffer[pixel_r, pixel_c] = depth
                image[pixel_r, pixel_c] = color
            if (d<=0):
                d += d_vertical
            else:
                d += d_diagonal
                pixel_c += step_col
                
##
# Given an array of points, a color and a type of geometry,
# draws that geometry by using draw_line function
##
def draw_geometry(image, points_array, color, type):
    for i in range(0, len(points_array)-1, 1):
        draw_line(image, points_array[i], points_array[i + 1], color)

    if (type == GEOMETRY_REGION):
        draw_line(image, points_array[-1], points_array[0], color)

##
# Returns a cube centralized in (0,0,0). The cube is denoted 
# by an array with 6 faces, each is an array of points
#   parameters[0] - stands for the size of the sides of the cube
#   parameters[1] - stands for whether the faces of the cube will be stylished with a X or not
##
def get_cube(parameters):

    l = float(parameters[0]) / 2
    is_fancy_cube = parameters[1] == '1'

    # Array for a simple cube (with no X). Explaining this array:
    # - General structure: [ face[0], face[1], face[2], face[3], face[4], face[05] ]
    # - face[i] structure: [ point[0], point[1], point[2], point[3] ]
    # - point[i] structure: [ x, y, z, 1.0 ]
    faces = np.array([[[-l,-l,-l,1.0],[l,-l,-l,1.0],[l,l,-l,1.0],[-l,l,-l,1.0]],
                    [[-l,-l,l,1.0],[l,-l,l,1.0],[l,l,l,1.0],[-l,l,l,1.0]],
                    [[-l,-l,-l,1.0],[-l,l,-l,1.0],[-l,l,l,1.0],[-l,-l,l,1.0]],
                    [[l,-l,-l,1.0],[l,l,-l,1.0],[l,l,l,1.0],[l,-l,l,1.0]],
                    [[-l,-l,-l,1.0],[-l,-l,l,1.0],[l,-l,l,1.0],[l,-l,-l,1.0]],
                    [[-l,l,-l,1.0],[-l,l,l,1.0],[l,l,l,1.0],[l,l,-l,1.0]]])

    # For stylish cubes (with an X on each face)
    if is_fancy_cube:
        for i in range (0, 6, 1):
            # Reordering the point list (p[1] and p[2])
            aux = np.copy(faces[i][1])
            faces[i][1] = np.copy(faces[i][2])
            faces[i][2] = np.copy(aux)
    
    return np.array(faces)

##
# Given an array of faces and a color
# draws a cube by using draw_geometry set to GEOMETRY_REGION
##
def draw_cube(image, faces, color):
    for i in range(0, len(faces), 1):
        draw_geometry(image, faces[i], color, GEOMETRY_REGION)

##
# Given an array of points and a color
# draws a sphere using draw_line
##
def draw_SPH(image, points_array, color):
    
    #desenha os meridianos
    for j in range(0,len(points_array)):  
        for i in range(0,len(points_array[j])-1): 
             draw_line(image, points_array[j][i], points_array[j][i+1],color)
        draw_line(image, points_array[j][-1], points_array[j][0],color)
    
    #desenha os paralelos
    for i in range(1,len(points_array)+1):    
        for j in range(0,(len(points_array)-1)):
                draw_line(image, points_array[j][i], points_array[j+1][i],color)
                draw_line(image, points_array[j][-i], points_array[j+1][-i],color)
        draw_line(image, points_array[-1][i], points_array[0][-i],color)
        draw_line(image, points_array[-1][-i], points_array[0][i],color)
    
 
# ---------- Main routine ----------

# Parses and checks command-line arguments
if len(sys.argv)!=3:
    print("usage: python draw_2d_model.py <input.dat> <output.ppm>\n"
          "       interprets the drawing instructions in the input file and renders\n"
          "       the output in the NETPBM PPM format into output.ppm")
    sys.exit(1)

input_file_name  = sys.argv[1]
output_file_name = sys.argv[2]

# Reads input file and parses its header
with open(input_file_name, 'rt', encoding='utf-8') as input_file:
    input_lines = input_file.readlines()

if input_lines[0] != 'EA979V4\n':
    print(f'input file format not recognized!', file=sys.stderr)
    sys.exit(1)

dimensions = input_lines[1].split()
width = int(dimensions[0])
height = int(dimensions[1])

if width<=0 or width>MAX_SIZE or height<=0 or height>MAX_SIZE:
    print(f'input file has invalid image dimensions: must be >0 and <={MAX_SIZE}!', file=sys.stderr)
    sys.exit(1)

# Creates image
image = np.full((height, width, CHANNELS_N), fill_value=DEFAULT_BACKGROUND, dtype=IMAGE_DTYPE)
zbuffer = np.full((height, width,), fill_value=ZBUFFER_BACKGROUND, dtype=ZBUFFER_DTYPE)
color = DEFAULT_COLOR

#
# TODO: Inicialize as demais variaveis
#
transformation_matrix = np.eye(4)
viewport_matrix = np.eye(4)
matrix_stack = []
# Main loop - interprets and renders drawing commands
for line_n,line in enumerate(input_lines[2:], start=3):

    if len(line)>MAX_LINE_LEN:
        print(f'line {line_n}: line too long!', file=sys.stderr)
        sys.exit(1)

    if not line.strip():
        # Blank line - skips
        continue
    if line[0] == '#':
        # Comment line - skips
        continue

    tokens = line.strip().split()
    command = tokens[0]
    parameters = tokens[1:]
    def check_parameters(n):
        if len(parameters) != n:
            print(f'line {line_n}: command {command} expected {n} parameters but got {len(parameters)}!',
                  file=sys.stderr)
            sys.exit(1)


    if command == SET_BACKGROUND_COLOR:
        # Clears with new background color
        check_parameters(CHANNELS_N)
        background_color = np.array(parameters, dtype=IMAGE_DTYPE)
        image[...] = background_color
        zbuffer[...] = ZBUFFER_BACKGROUND
        
    elif command == SET_COLOR:
        ##
        # Changes the color used to draw lines
        ##
        check_parameters(CHANNELS_N)
        color = np.array(parameters, dtype=IMAGE_DTYPE)
        
    elif command == DRAW_GEOMETRY_LINE:
        ##
        # Draws a line based on 2 input points
        ##
        check_parameters(6)
        points_h= get_points(parameters)

        # Applies the matriz operator to the coordenates
        points_ht = transform(transformation_matrix,points_h)
        
        # Applies the perspective matriz to the coordenates
        points_hp = transform(viewport_matrix,points_ht)

        # Draws a line
        draw_line(image, points_hp[0], points_hp[1], color)
   
    elif command == DRAW_GEOMETRY_POLYLINE:
        ##
        # Draws a polyline based on N input points
        ##
        check_parameters(3*(int(parameters[0])) +1)
        points_h = get_points(parameters)
    
        # Applies the matrix operator to the coordinates
        points_ht = transform(transformation_matrix,points_h)
            
        # Applies the perspective matriz to the coordenates
        points_hp = transform(viewport_matrix,points_ht)
    
        # Draws a polyline
        draw_geometry(image, points_hp, color, GEOMETRY_POLYLINE)

    elif command == DRAW_GEOMETRY_REGION:
        ##
        # Draws a region based on N input points
        # A line will be drawn from the last point to the first one to enclose the region
        ##
        check_parameters(3*(int(parameters[0])) +1)
        points_h = get_points(parameters)

        # Applies the matrix operator to the coordinates
        points_ht = transform(transformation_matrix,points_h)
        
        # Applies the perspective matriz to the coordenates
        points_hp = transform(viewport_matrix,points_ht)

        # Draws a polygon
        draw_geometry(image, points_hp, color, GEOMETRY_REGION)
    
    
    elif command == SET_TRANSFORMATION_MATRIX:
        ##
        # Replaces the current transformation matrix applied to all input geometries
        ##
        check_parameters(16)

        # Replacing the current transformation matrix with the input
        transformation_matrix = get_matrix(parameters)    
        
    elif command == SET_VIEWPORT_MATRIX:
        ##
        # Replaces the current perspective matrix applied to all input geometries
        ##
        check_parameters(16)

        # Replacing the current perspective matrix with the input
        viewport_matrix = get_matrix(parameters)    
        
    elif command == APPLY_TRANSFORMATION_MATRIX:
        ##
        # Applies a new transformation matrix to the current transformation matrix
        ##
        check_parameters(16)

        # New transformation matrix is the current multiplied by the inputted one
        input_matrix = get_matrix(parameters)
        transformation_matrix = np.matmul(transformation_matrix,input_matrix)

    elif command == DRAW_GEOMETRY_SPHERE :
        ##
        # Draws a sphere according to inputted radius, and number of meridians and parallels
        ##

        check_parameters(3)
        
        points_SPH = get_sph(float(parameters[0]),int(parameters[1]),int(parameters[2]))
        
        # Applies the matrix operator to the coordinates
        points_SPHt = transform_sph(transformation_matrix,points_SPH)
            
        # Applies the perspective matrix to the coordenates
        points_SPHp = transform_sph(viewport_matrix,points_SPHt)
        
        draw_SPH(image, points_SPHp, color)

    elif command == DRAW_GEOMETRY_CUBE:
        ##
        # Draws a cube according to inputted size and style (1 for crossed faces)
        ##

        check_parameters(2)

        # Gets the 6 faces to be drawn
        faces = get_cube(parameters)

        # Applies transformations in every face
        for i in range(0, len(faces), 1):
            faces[i] = transform(viewport_matrix,transform(transformation_matrix,faces[i]))
            
        # Draws the cube
        draw_cube(image, faces, color)
        
    elif command == PUSH_TRANSFORMATION_MATRIX:
        ##
        # Pushes current transoformation matrix to a stack
        ##
        matrix_stack.append(transformation_matrix)

    elif command == POP_TRANSFORMATION_MATRIX:
        ##
        # Replaces the current transformation matrix to a new one popped from a stack
        ##
        transformation_matrix = matrix_stack.pop()

    else:
        print(f'line {line_n}: unrecognized command "{command}"!', file=sys.stderr)
        sys.exit(1)

# If we reached this point, everything went well - outputs rendered image file
with open(output_file_name, 'wb') as output_file:
    save_ppm(image, output_file)
