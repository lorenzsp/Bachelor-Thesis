%% binary neutron star
video = VideoWriter('bns.mp4','MPEG-4'); %create the video object
open(video); %open the file for writing
N=335;
for ii=1:N %where N is the number of images
  filename = sprintf('rho_000000%03d.png',ii);
  I = imread(filename); %read the next image
  writeVideo(video,I); %write the image to file
end
close(video); %close the file