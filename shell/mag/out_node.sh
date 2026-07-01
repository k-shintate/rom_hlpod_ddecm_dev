
#!/bin/bash

for i in $(seq 0 8); do
  echo "processing partition ${i}"

  python3 ../../shell/mag/out_node.py \
    --elem parted.0/nedelec_elem.dat.${i} \
    --coord parted.0/node_coordinate_elem.dat.${i} \
    --output parted.0/sorted_nodes.dat.${i} \
    --output-local-elem parted.0/sorted_local_elem.dat.${i} \
    --map-output parted.0/global_local_map.dat.${i}
done
