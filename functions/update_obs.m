function [obs] = update_obs(obs,dt)
num_obs = size(obs,2);
for j = 1:num_obs
    switch obs{j}.shape
        case 'circle'
        obs{j}.c = obs{j}.c+obs{j}.v*dt;
        case 'rectangle'
        obs{j}.c = obs{j}.c+obs{j}.v*dt;        
        obs{j}.poly = [obs{j}.c(1)-obs{j}.w/2  obs{j}.c(1)+obs{j}.w/2  obs{j}.c(1)+obs{j}.w/2  obs{j}.c(1)-obs{j}.w/2;
                       obs{j}.c(2)-obs{j}.l/2  obs{j}.c(2)-obs{j}.l/2  obs{j}.c(2)+obs{j}.l/2  obs{j}.c(2)+obs{j}.l/2];
    end
end

end
