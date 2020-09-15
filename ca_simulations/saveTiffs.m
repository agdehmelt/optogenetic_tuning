if(save_tiffs==true)
   if save_1==true;
       div=single(squeeze(O(:,:,1,1,save_start:10:save_end)));
       tif32write(div,[output_folder 'Ga_' num2str(sim_num) '.tif']);
       fprintf(['saved tif stack: ' output_folder 'Ga_' num2str(sim_num) '.tif']);
   end
   if save_2==true;
       div=single(squeeze(O(:,:,1,2,save_start:10:save_end)));
       tif32write(div,[output_folder 'Gi_' num2str(sim_num) '.tif']);
       fprintf(['saved tif stack: ' output_folder 'Gi_' num2str(sim_num) '.tif']);
   end
   if save_3==true;
       div=single(squeeze(O(:,:,1,3,save_start:10:save_end)));
       tif32write(div,[output_folder 'Ma_' num2str(sim_num) '.tif']);
       fprintf(['saved tif stack: ' output_folder 'Ma_' num2str(sim_num) '.tif']);
   end
   if save_4==true;
       div=single(squeeze(O(:,:,1,4,save_start:10:save_end)));
       tif32write(div,[output_folder 'Mi_' num2str(sim_num) '.tif']);
       fprintf(['saved tif stack: ' output_folder 'Mi_' num2str(sim_num) '.tif']);
   end
   if save_5==true;
       div=single(squeeze(O(:,:,1,5,save_start:10:save_end)));
       tif32write(div,[output_folder 'Ra_' num2str(sim_num) '.tif']);
       fprintf(['saved tif stack: ' output_folder 'Ra_' num2str(sim_num) '.tif']);
   end
   if save_6==true;
       div=single(squeeze(O(:,:,1,6,save_start:10:save_end)));
       tif32write(div,[output_folder 'Ri_' num2str(sim_num) '.tif']);
       fprintf(['saved tif stack: ' output_folder 'Ri_' num2str(sim_num) '.tif']);
   end
end