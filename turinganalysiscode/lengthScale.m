filename = 'New Data/crop3.bmp';
image = imread(filename);
figure, imshow(image), title 'original';
bw = im2bw(image, graythresh(image)+.15);

[x, y]=size(bw);
halfx = x/2;
halfy = y/2;

bw = bwareaopen(bw, 50);
figure, imshow(bw), title 'black and white';
cc=bwconncomp(bw);
stats = regionprops(cc, 'Centroid','Area','MinorAxisLength','MajorAxisLength');
spotSize = [stats.Area];
minor=[stats.MinorAxisLength];
major=[stats.MajorAxisLength];
sizeArray=length(stats);
Distance = zeros(1, sizeArray);
minDistance = zeros(1, sizeArray);
for i=1:length(stats)
    for j=1:length(stats)
        a=stats(i).Centroid-stats(j).Centroid;
        diffx = abs(a(1));
        diffy = abs(a(2));
        if diffx > halfx
            diffx = x-diffx;
        end
        if diffy > halfy
            diffy = y - diffy;
        end
        dist=sqrt((diffx)^2+(diffy)^2);
        if(dist==0)
            dist=1000;
        end
        Distance(j)=dist;
    end
    minDistance(i)=min(Distance);
end
figure, hist(minDistance,40), title 'Histogram of minimium distance to next spot', xlabel 'distance to next spot (px)', ylabel 'frequency';
a=findobj(gcf); % get the handles associated with the current figure

allaxes=findall(a,'Type','axes');
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');

set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,...
    'FontSize',14);
set(alllines,'Linewidth',2);
set(alltext,'FontName','Arial','FontWeight','Bold','FontSize',14);
figure, histfit(minDistance), title 'best fit normal distribution';
pd=fitdist(minDistance(:),'normal')
mean(minDistance)
std(minDistance)
median(minDistance)
paramEstsMinima = evfit(-minDistance)
y = linspace(15,80,65);
figure, hist(minDistance,15), title 'best fit Gumbel';
p = evpdf(-y,paramEstsMinima(1),paramEstsMinima(2));
line(y,3*length(minDistance)*p,'color','r');
sizeSpot=sqrt(spotSize/pi);
figure, histfit(sizeSpot), title 'best fit normal size of spot'
pd2=fitdist(sizeSpot','normal')
figure, hist(sizeSpot,40), title 'histogram of sqrt area';
figure, hist(minor, 40), title 'minor axis length'
figure, hist(major, 40), title 'major axis length'