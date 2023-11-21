# Mesh Clustering, Reader and Writer

As per the first step in the Nanite pipeline, this project splits up a mesh into clusters. The algorithm attempts to keep the clusters as localized as possible (minimizing strips) in order to keep the rendering efficient. This is then run recursively and the mesh encoded to create nanite like meshes that I use in my own recreation of the technology.

## Example Output on a Bunny Mesh in Blender

![image](https://github.com/mwjrink/gim/assets/29183162/23ad4e2c-4f85-4347-9621-c2e129f3d509)
