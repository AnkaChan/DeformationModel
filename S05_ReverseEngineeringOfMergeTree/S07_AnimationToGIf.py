import imageio, glob

if __name__ == '__main__':
    filenames = glob.glob('./Data/S06_TransitionAnimartion.py/Animation/*.png')
    images = []
    for filename in filenames:
        images.append(imageio.imread(filename))
    imageio.mimsave('./Data/S06_TransitionAnimartion.py/movie.gif', images)