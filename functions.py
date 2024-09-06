import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import imageio
matplotlib.use('TkAgg')

def read_outputdt(folder_path):
    file_name = os.path.join(folder_path, 'runscript.txt')
    df = np.genfromtxt(file_name)
    outputline = df[6]
    output_dt = outputline[2]
    return output_dt

def write_runscript(lines, folder_path):
    # lines = ['Readme', 'How to write text files in Python']
    file_name = os.path.join(folder_path, 'runscript.txt')
    with open(file_name, 'w') as f:
        for line in lines:
            f.write(line)
            f.write('\n')

def save_frame(mat_list, total, folder_path, output_step):
    # <--------------------------------------------------------------------------------------------------------
    # Density[Total] folder
    frame_folder = os.path.join(folder_path, 'Density[Total]')
    if not os.path.exists(frame_folder):
        os.mkdir(frame_folder)
    file_path = os.path.join(frame_folder, 'Density[Total]_{}.png'.format(output_step))
    fig, ax = plt.subplots()
    im = ax.imshow(total.density)
    fig.colorbar(im, label="Density", orientation="vertical")
    im.set_clim(0,mat_list[0].material.density0)
    fig.savefig(file_path, bbox_inches='tight')
    plt.close(fig)
    # <--------------------------------------------------------------------------------------------------------
    # Momentumx[Total] folder
    frame_folder = os.path.join(folder_path, 'Momentumx[Total]')
    if not os.path.exists(frame_folder):
        os.mkdir(frame_folder)
    file_path = os.path.join(frame_folder, 'Momentumx[Total]_{}.png'.format(output_step))
    fig, ax = plt.subplots()
    im = ax.imshow(total.momentum_x, cmap='jet')
    fig.colorbar(im, label="Momentumx", orientation="vertical")
    im.set_clim(-1.0e+6, 0.0)
    fig.savefig(file_path, bbox_inches='tight')
    plt.close(fig)
    # <--------------------------------------------------------------------------------------------------------

def make_gifs(folder_path):
    # <--------------------------------------------------------------------------------------------------------
    frame_folder = os.path.join(folder_path, 'Density[Total]')
    contents = os.listdir(frame_folder)
    image_paths = []
    for image in contents:
        image_paths.append(os.path.join(frame_folder, image))
    images = []
    image_paths = sorted(image_paths, key=lambda y: int(y.split("Density[Total]_")[1].split(".png")[0]))
    for filename in image_paths:
        images.append(imageio.v2.imread(filename))
    gif_path = os.path.join(folder_path, 'Density[Total].gif')
    imageio.mimsave(gif_path, images, duration=0.5)
    # <--------------------------------------------------------------------------------------------------------
    frame_folder = os.path.join(folder_path, 'Momentumx[Total]')
    contents = os.listdir(frame_folder)
    image_paths = []
    for image in contents:
        image_paths.append(os.path.join(frame_folder, image))
    images = []
    image_paths = sorted(image_paths, key=lambda y: int(y.split("Momentumx[Total]_")[1].split(".png")[0]))
    for filename in image_paths:
        images.append(imageio.v2.imread(filename))
    gif_path = os.path.join(folder_path, 'Momentumx[Total].gif')
    imageio.mimsave(gif_path, images, duration=0.5)
    # <--------------------------------------------------------------------------------------------------------
    print('Animations saved')