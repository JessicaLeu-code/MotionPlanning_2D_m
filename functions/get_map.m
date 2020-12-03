function [grid] =get_map(obs, range)
grid_x = linspace(range(1),range(2),round((range(2)-range(1))/0.001));
grid_y = linspace(range(3),range(4),round((range(4)-range(3))/0.001));
grid= false(size(grid_y,2),size(grid_x,2));

for i= 1:size(grid_x,2)
    for y = 1:size(grid_y,2)        
        for j = 1:size(obs,2)
               if ~isempty(obs{j})
                   switch obs{j}.shape
                       case 'circle' 
                           if i==300
                               eee= 1;
                           end
                           if norm([grid_x(i);grid_y(y)] - obs{j}.c)<obs{j}.m
                              grid(size(grid_y,2)-y+1,i) = true;
                              break;                           
                           end          
                       case 'rectangle'
                           region_obs = [obs{j}.w/2;obs{j}.l/2];
                           if (obs{j}.c-region_obs)<[grid_x(i);grid_y(y)]
                               if [grid_x(i);grid_y(y)]<(obs{j}.c+region_obs)
                                   grid(size(grid_y,2)-y+1,i) = true;
                                   break;
                               
                               end                           
                           end 
                    end
                   
               end
         end
    end
end

        
