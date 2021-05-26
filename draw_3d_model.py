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
SET_BACKGROUND_COLOR = 'c'
SET_COLOR = 'C'
INPUT_GEOMETRY_LINE = 'L'
INPUT_GEOMETRY_POLYLINE = 'P'
INPUT_GEOMETRY_REGION = 'R'
SET_TRANSFORMATION_MATRIX = 'M'
APPLY_TRANSFORMATION_MATRIX = 'm'


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
    # separates the input form of P comand, and R and P comand
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



# ---------- Drawing/model routines ----------
def draw_line(image, V0, V1, color):
    #Transform the vectors in homogeneous coordinates to R3 int
    x0 = round((V0[0]/V0[3])+width/2)
    x1 = round((V1[0]/V1[3])+width/2)
    y0 = round((V0[1]/V0[3])+height/2)
    y1 = round((V1[1]/V1[3])+height/2)
    # Adicionar Zbuffer
    z0 = round(V0[2])
    z1 = round(V1[2])
    #assert x0 >=0 and x0 < width and x1 >=0 and x1 < width and y0 >=0 and y0 < height and y1 >=0 and y1 < height
    # Computes differences
    dx = x1-x0
    dy = y1-y0
    dc = abs(dx) # delta x in book - here we are using row, col coordinates
    dr = abs(dy) # delta y in book
    if dr <= dc:
        # Line inclination is at most 1
        # Swaps points if c1<c0 and converts x,y coordinates to row,col coordinates
        # dx>=0 => x1>=x0 => c1>=x0
        r0 = height-1-y0 if dx>=0 else height-1-y1
        r1 = height-1-y1 if dx>=0 else height-1-y0
        c0 =     x0 if dx>=0 else x1
        c1 =     x1 if dx>=0 else x0
        # Implements Bresenham's midpoint algorithm for lines
        # (Klawonn. Introduction to Computer Graphics. 2nd Edition. Section 4.2, pp. 45–53)
        # ...deltas of Bressenham's algorithm
        d_horizontal = 2*dr      # delta east in book
        d_diagonal   = 2*(dr-dc) # delta northeast in book
        # ...draws line
        pixel_r = r0
        step_row = 1 if r1>=r0 else -1
        d = 2*dr - dc # starting D value, D_init in book
        for pixel_c in range(c0, c1+1):
            depth = z0 + (z1-z0)*(pixel_c-c0)/(c1-c0)
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
        # (Klawonn. Introduction to Computer Graphics. 2nd Edition. Section 4.2, pp. 45–53)
        # ...deltas of Bressenham's algorithm - same as above, but with coordinates inverted
        d_vertical = 2*dc
        d_diagonal = 2*(dc-dr)
        pixel_r = r0
        pixel_c = c0
        step_col = 1 if c1>=c0 else -1
        d = 2*dc - dr # starting D value, D_init in book
        for pixel_r in range(r0, r1+1):
            depth = z0 + (z1-z0)*(pixel_r-r0)/(r1-r0)
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
    


# ---------- Main routine ----------

# Parses and checks command-line arguments
if len(sys.argv)!=3:
    print("usage: python draw_2d_model.py <input.dat> <output.ppm>\n"
          "       interprets the drawing instructions in the input file and renders\n"
          "       the output in the NETPBM PPM format into output.ppm")
    sys.exit(1)

input_file_name  = sys.argv[1]
output_file_name = sys.argv[2]
matrix_stack = []

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


    if command == 'c':
        # Clears with new background color
        check_parameters(CHANNELS_N)
        background_color = np.array(parameters, dtype=IMAGE_DTYPE)
        image[...] = background_color
        zbuffer[...] = ZBUFFER_BACKGROUND
        
    elif command == "C":
        ##
        # Changes the color used to draw lines
        ##
        check_parameters(CHANNELS_N)
        color = np.array(parameters, dtype=IMAGE_DTYPE)
        
    elif command == "L":
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
   
    elif command == "P":
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

    elif command == "R":
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
    
    
    elif command == "M":
        ##
        # Replaces the current transformation matrix applied to all input geometries
        ##
        check_parameters(16)

        # Replacing the current transformation matrix with the input
        transformation_matrix = get_matrix(parameters)    
        
    elif command == "V":
        ##
        # Replaces the current perspective matrix applied to all input geometries
        ##
        check_parameters(16)

        # Replacing the current perspective matrix with the input
        viewport_matrix = get_matrix(parameters)    
        
    elif command == "m":
        ##
        # Applies a new transformation matrix to the current transformation matrix
        ##
        check_parameters(16)

        # New transformation matrix is the current multiplied by the inputted one
        input_matrix = get_matrix(parameters)
        transformation_matrix = np.matmul(transformation_matrix,input_matrix)
      
    elif command == "PUSH":
        matrix_stack.append(transformation_matrix)
        print("Comando push")

    elif command == "POP":
        print("comando POP")
        transformation_matrix = matrix_stack.pop()

    elif command == "SPH":
        print("comando SPH")

    elif command == "CUB":
        print("comando CUB")

        
    #
    # TODO: Implemente os demais comandos
    #

    else:
        print(f'line {line_n}: unrecognized command "{command}"!', file=sys.stderr)
        sys.exit(1)

# If we reached this point, everything went well - outputs rendered image file
with open(output_file_name, 'wb') as output_file:
    save_ppm(image, output_file)
