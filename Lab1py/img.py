from PIL import Image
import numpy as np
from math import cos, sin, sqrt

def dotted_line(image, x0, y0, x1, y1, color):
    count = sqrt((x0 - x1)**2 + (y0 - y1)**2)
    step = 1.0/count
    for t in np.arange (0, 1, step):
        x = round ((1.0 - t)*x0 + t*x1)
        y = round ((1.0 - t)*y0 + t*y1)
        image[y, x] = color

def draw_line(image, x0, y0, x1, y1, color):
    count = 100
    step = 1.0/count
    for t in np.arange (0, 1, step):
        x = round ((1.0 - t)*x0 + t*x1)
        y = round ((1.0 - t)*y0 + t*y1)
        image[y, x] = color

def x_loop_line(image, x0, y0, x1, y1, color):
    xchange = False
    
    if (abs(x0 - x1) < abs(y0 - y1)):
        x0, y0 = y0, x0
        x1, y1 = y1, x1
        xchange = True
    if (x0 > x1):
        x0, x1 = x1, x0
        y0, y1 = y1, y0

    y = y0
    dy = 2*abs(y1 - y0)
    derror = 0.0
    y_update = 1 if y1 > y0 else -1
            
    for x in range (x0, x1):
        t = (x - x0)/(x1 - x0)
        y = round((1.0 - t)*y0 + t*y1)
        
        if (xchange):
            image[x, y] = color
        else:
            image[y, x] = color

        derror += dy
        if (derror > (x1 - x0)):
            derror -= (x1 - x0)
            y += y_update


img1 = np.zeros((200, 200, 3), dtype = np.uint8)
img2 = np.zeros((200, 200, 3), dtype = np.uint8)
img3 = np.zeros((200, 200, 3), dtype = np.uint8)

for k in range(13):
    x0,y0 = 100,100
    x1= 100+95*cos((2*np.pi/13)*k)
    y1= 100+95*sin((2*np.pi/13)*k)
    dotted_line(img1, x0, y0, x1, y1, 255)
    draw_line(img2, x0, y0, x1, y1, 255)
    x_loop_line(img3, x0, y0, int(x1), int(y1), 255)

#for i in range (600):
 #   for j in range (800):
  #      if (i <j):   
   #         img[i,j] = [0, 200, 200]
Image.fromarray(img1).save("img1.png")
Image.fromarray(img2).save("img2.png")
Image.fromarray(img3).save("img3.png")
