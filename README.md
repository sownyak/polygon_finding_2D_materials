# polygon_finding_2D_materials
#You can calculate the polygons of your 2D materials (periodic) structres :
# This Code Can be used for Graphene, Graphene Oxide, B-N doped or any other 2D materials.
# Use the bond information of your structure (LAMMPS generate automatically) and generate an file with the image information of your structure. then sort the image file. 
#run the following command in your terminal. "cycle_file.txt" will be your output file along with a coordinate file.

python3 polygons.py image_final.txt structure.xyz cycle_file.txt



# Cite our paper :
@article{mondal2024quantifying,
  title={Quantifying defects in graphene oxide structures},
  author={Mondal, Sownyak and Ghosh, Soumya},
  journal={Carbon Trends},
  volume={14},
  pages={100323},
  year={2024},
  publisher={Elsevier}
}