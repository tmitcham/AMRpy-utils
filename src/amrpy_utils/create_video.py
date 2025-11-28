import glob
from pathlib import Path
import shutil
import cv2


def create_video(output_path, input_patterns, fps):

    # Create a list of all the input image files
    FILES = []
    for i in input_patterns:
        FILES += glob.glob(i)

    # Get the filename from the output path
    filename = Path(output_path).name
    print(f'Creating video "{filename}" from images "{FILES}"')

    # Load the first image to get the frame size
    frame = cv2.imread(FILES[0])
    height, width, layers = frame.shape

    # Create a VideoWriter object to write the video file
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    video = cv2.VideoWriter(filename=filename, fourcc=fourcc, fps=fps, frameSize=(width, height))

    # Loop through the input images and add them to the video
    for image_path in FILES:
        print(f'Adding image "{image_path}" to video "{output_path}"... ')
        video.write(cv2.imread(image_path))

    # Release the VideoWriter and move the output file to the specified location
    cv2.destroyAllWindows()
    video.release()

    shutil.move(filename, output_path)