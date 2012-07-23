##############################################################################################################
# -*- coding: cp1252 -*-
# 
# This file is part of the mRF experiment
# Date: May 2011
# Author: Franck DERNONCOURT <franck.dernoncourt@gmail.com>
# License: GPLv3
# ----
# Python: 2.6.6 32bits
# Pygame: 1.9.1
#
##############################################################################################################


import os, sys, pygame, random
from pygame.locals import *

# Global constants
window_size = (500, 600)
window_center = (window_size[0] / 2.0, window_size[1] / 2.0)
window_scale = min(window_size) / 2.0
black_color = (0, 0, 0)
grey_color = (128, 128, 128)
white_color = (255, 255, 255)
robot_color = (0, 0, 255)
robot_radius = 18
map_size = (400, 400)
map_center_location = ((window_size[0] - map_size[0]) / 2, 20)
nb_tiles_h = 5;
nb_tiles_w = 5;
tile_size_real_h = map_size[1] / nb_tiles_h
tile_size_real_w = map_size[0] / nb_tiles_h;
tile_size_pixel_h = tile_size_real_h
tile_size_pixel_w = tile_size_real_w
fps = 30
text_y = map_center_location[1] + map_size[1]
text_h = 20 


def coord(real):
    ''' 
    Takes real coordinates, returns pixel coordinates
    '''
    return (int(round(real[0] * window_scale + window_center[0])), \
            int(round(-real[1] * window_scale + window_center[1])))


def coord_map(coord):
    ''' 
    Takes map coordinates, returns real coordinates
    '''
    return (coord[0] + map_center_location[0], coord[1] + map_center_location[1])
    
def draw_tiles(surf):
    ''' 
    Draw tiles (black and white)
    WARNING: make sure the tiles are the same than in mrf.cpp's tiles_pos
    variable.
    '''
    tiles_pos = dict()
    tiles_pos[(1, 2)] = 'white'
    tiles_pos[(3, 3)] = 'white'
    tiles_pos[(3, 1)] = 'black'
    tiles_pos[(2, 3)] = 'black'
    for k, v in tiles_pos.items():
        if (v == 'white'):
            rect_color = white_color
        elif (v == 'black'):
            rect_color = black_color

        pygame.draw.rect(surf, rect_color, \
                        pygame.Rect(coord_map((tile_size_pixel_w * k[1], tile_size_pixel_h * k[0])),\
                                    (tile_size_pixel_w, tile_size_pixel_h)))

def initialize_map(surf, fontobject):
    ''' 
    Draw tiles (black and white)
    WARNING: make sure the tiles are the same than in mrf.cpp's tiles_pos
    variable.
    '''    
    # Draw background
    surf.fill(black_color)
    
    # Draw grid
    pygame.draw.rect(surf, grey_color, pygame.Rect(map_center_location, map_size))
    draw_tiles(surf)
    
    # Draw variable names
    surf.blit(fontobject.render('action:', 1, (255, 255, 0)),
                (3, text_y))
    
    surf.blit(fontobject.render('e =', 1, (255, 255, 0)),
                (27, text_y + text_h))
    
    surf.blit(fontobject.render('pe =', 1, (255, 255, 0)),
                (20, text_y + text_h * 2))
    
    surf.blit(fontobject.render('bl =', 1, (255, 255, 0)),
                (20, text_y + text_h * 3))
            
    surf.blit(fontobject.render('br =', 1, (255, 255, 0)),
                (20, text_y + text_h * 4))
            
    surf.blit(fontobject.render('ld =', 1, (255, 255, 0)),
                (20, text_y + text_h * 5))
            
    surf.blit(fontobject.render('lb =', 1, (255, 255, 0)),
                (20, text_y + text_h * 6))

def explicit_action(selected_action):
    if (selected_action == 0):
        return 'wander'
    elif (selected_action == 1):
        return ' avoid'
    elif (selected_action == 2):
        return ' reload_on_dark'
    elif (selected_action == 3):
        return ' reload_on_light'
    elif (selected_action == 4):
        return ' rest'
    
    return 'WARNING: No action detected'
    
def main():
    '''
    This is the main function
    '''
    global fps
    try:

        # Graphics initialization
        full_screen = False
        pygame.init()
        if full_screen:
            surf = pygame.display.set_mode(window_size, HWSURFACE | FULLSCREEN | DOUBLEBUF)
        else:
            surf = pygame.display.set_mode(window_size)
            
        # Initialize variables text
        fontobject = pygame.font.SysFont('Arial', 18)

        pygame.display.set_caption("mRF robot simulation replay")
        initialize_map(surf, fontobject)
        
        # Draw robot
        pygame.draw.circle(surf, robot_color, coord_map((200, 200)), robot_radius, 0)
        
        # Draw step number 
        surf.blit(fontobject.render('step#: ' + str(0), 1, (255, 255, 0)),
                (map_center_location[0] + map_size[0] - 100, text_y))
    
        # Update display
        pygame.display.flip()
        
        # Open log file
        moves_file = open('moves3.txt', 'r')
        
        # Wait for keydown before beginning the video
        no_keydown = True
        while no_keydown:
            for event in pygame.event.get():
                if (event.type == KEYDOWN):
                    no_keydown = False

        # Main loop
        num_step = 0
        while True:
            # Read line
            moves = moves_file.readline().split(',')
            num_step = float(moves[0])
            robot_x = float(moves[1])
            robot_y = float(moves[2])
            selected_action = float(moves[3])
            e = float(moves[4])
            pe = float(moves[5])
            bl = float(moves[6])
            br = float(moves[7])
            ld = float(moves[8])
            lb = float(moves[9])
            
            # Draw robot
            initialize_map(surf, fontobject)
            pygame.draw.circle(surf, robot_color, coord_map((int(robot_y), int(robot_x))), robot_radius, 0)
            
            # Draw variables text            
            surf.blit(fontobject.render(explicit_action(selected_action), 1, (255, 255, 0)),
                (50, text_y))
            
            surf.blit(fontobject.render(str(e), 1, (255, 255, 0)),
                (50, text_y + text_h))
            
            surf.blit(fontobject.render(str(pe), 1, (255, 255, 0)),
                (50, text_y + text_h * 2))
            
            surf.blit(fontobject.render(str(bl), 1, (255, 255, 0)),
                (50, text_y + text_h * 3))
            
            surf.blit(fontobject.render(str(br), 1, (255, 255, 0)),
                (50, text_y + text_h * 4))
            
            surf.blit(fontobject.render(str(ld), 1, (255, 255, 0)),
                (50, text_y + text_h * 5))
            
            surf.blit(fontobject.render(str(lb), 1, (255, 255, 0)),
                (50, text_y + text_h * 6))
            
            # Draw step number 
            surf.blit(fontobject.render('step#: ' + str(num_step), 1, (255, 255, 0)),
                (map_center_location[0] + map_size[0] - 100, text_y))
            
            # Draw saliences
            sep = ' '
            surf.blit(fontobject.render('saliences: ' + str(moves[10]) + sep + str(moves[11]) + sep + str(moves[12]) + sep + str(moves[13]) , 1, (255, 255, 0)),
                (map_center_location[0] + map_size[0] - 300, text_y + 20))
            
            # Draw output
            surf.blit(fontobject.render('outputs: ' + str(moves[15]) + sep + str(moves[16]) + sep + str(moves[17]) + sep + str(moves[18]) , 1, (255, 255, 0)),
                (map_center_location[0] + map_size[0] - 300, text_y + 40))
            
            # Draw contrast
            surf.blit(fontobject.render('contrast: ' + str(moves[14]), 1, (255, 255, 0)),
                (map_center_location[0] + map_size[0] - 300, text_y + 60))
            
            # Draw FPS
            surf.blit(fontobject.render('FPS: ' + str(fps), 1, (255, 255, 0)),
                (window_size[0] - 70, window_size[1] - 20))
            
            # Handle events during experience            
            for event in pygame.event.get():
                if (event.type == KEYDOWN and event.key == K_ESCAPE)\
                    or event.type in (QUIT, MOUSEBUTTONDOWN):
                    pygame.quit()
                    sys.exit()
                if (event.type == KEYDOWN and (event.key == K_i)):
                    fps = fps + 5
                if (event.type == KEYDOWN and (event.key == K_d)):
                    fps = fps - 5
                    if (fps < 1 ):
                        fps = 1

            pygame.display.flip()
            pygame.time.delay(int(1.0/fps*1000))
            #num_step += 1

    finally:
        pygame.quit()
        #result_file.close()
        

if __name__ == "__main__":
    main()