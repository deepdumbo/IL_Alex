%{
    last update 2020/03/06, Hossein

    a function that shows a bunch of slides within an QSM image next to one another
%}

%function[]=Looped_Plotter(images)

% This function will display all of your images side by side, as the
% picking parts looped code does

%load your maps cell that you want to plot


set(gcf, 'Position', get(0,'Screensize')) %Maximize the image 
% figure('Position', [0,0,1920,1080])

num_loops=size(maps,1);

for i = 1:num_loops
    subplot (3,4,i)
    imagesc(maps{i,1}(:,:,80),[-0.3,0.3])
    axis image
    title(sprintf('TKD QSM, threshold %0.3f',maps{i,2}))
    colorbar
    impixelinfo
    colormap(bone)
end






% saveas(gcf, ['./Plots/',num2str(images.MaskErode) '_cube_' images.UnwrapType '_' images.RDFtype '_' images.QSMtype '_' images.orient])
%print(gcf, ['./Plots/',num2str(images.MaskErode) '_cube_' images.UnwrapType '_' images.RDFtype '_' images.QSMtype '_' images.orient],'-dpng','-r600')

%     parula     - Blue-green-orange-yellow color map
%     hsv        - Hue-saturation-value color map.x`
%     hot        - Black-red-yellow-white color map.
%     gray       - Linear gray-scale color map.
%     bone       - Gray-scale with tinge of blue color map.
%     copper     - Linear copper-tone color map.
%     pink       - Pastel shades of pink color map.
%     white      - All white color map.
%     flag       - Alternating red, white, blue, and black color map.
%     lines      - Color map with the line colors.
%     colorcube  - Enhanced color-cube color map.
%     vga        - Windows colormap for 16 colors.
%     jet        - Variant of HSV.
%     prism      - Prism color map.
%     cool       - Shades of cyan and magenta color map.
%     autumn     - Shades of red and yellow color map.
%     spring     - Shades of magenta and yellow color map.
%     winter     - Shades of blue and green color map.
%     summer     - Shades of green and yellow color map.
