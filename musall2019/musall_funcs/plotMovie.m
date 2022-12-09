function plotMovie(mov)
% mov (x,y,time)

% temp = reshape(mov, [], size(mov,3));
% % temp = double(temp) - (mean(temp,2));
% temp = temp - uint8(mean(temp,2));
% temp = reshape(temp,size(mov,1),size(mov,2),size(mov,3));

temp = mov;

figure; 
for i = 1:size(mov,3)
    imagesc(temp(:,:,i))
    colormap(gray)
    drawnow
%     pause(0.05)
end

end