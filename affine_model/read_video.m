
%% reading the videos
dir_video_input='input_video';
videos=dir(dir_video_input);
disp(videos)

%% processing each video 

for i=1:length(videos)
  video_name = videos(i).name;
  angle=sscanf(video_name,'angle_%d');
  disp(videos(i).name)
  
  if(~isempty(angle))
        path=strcat('input_video/',video_name);
        v_obj1 = VideoReader(path);
        v_obj = VideoReader(path); %copy for getting the number of frames
        h = v_obj.Height;
        w = v_obj.Width;
        f_rate = v_obj.FrameRate;
        clc;
        
        %% select the input features
        grand_angle=input('Type the value of the panoramic angle   ');
        angle_sample=360;
        
        while(angle_sample>grand_angle/2)
        angle_sample=input('Type the value of the distance-angle from one frame to the other   ');
        end
        clc;
        
        %% reading the video
        rate=ceil(grand_angle/angle_sample);
        n_frames = v_obj1.NumberOfFrames;

        %% subsampling the video to extract the frames
        space_frame=floor(n_frames/rate);
        cont=1;
        angle_str=num2str(angle);
        dir_name=strcat('input/angle_',angle_str);
        mkdir(dir_name);
        for j=1:n_frames
           if( rem(j,space_frame)==0)
               mov_final(:,:,:,cont)=readFrame(v_obj);
               
               if(cont<10) 
               cont_str=num2str(cont);
               cont_str=strcat('0',cont_str);
               end
               pathname=strcat(dir_name,'/i',cont_str,'.jpg');
               %disp(strcat('sampled frame number=',num2str(j)));
               imwrite( mov_final(:,:,:,cont),pathname);
               cont=cont+1;
           else
               readFrame(v_obj); 
           end
        end
  end
 
end
