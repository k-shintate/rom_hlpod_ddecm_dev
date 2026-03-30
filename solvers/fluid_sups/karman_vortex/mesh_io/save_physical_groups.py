import gmsh
import os
import logging
import sys

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

if len(sys.argv) < 2:
    logging.error("Please provide a mesh file as a command-line argument.")
    sys.exit(1)

mesh_file = sys.argv[1]

# Check if the file exists
if not os.path.isfile(mesh_file):
    logging.error(f"Mesh file '{mesh_file}' not found.")
    sys.exit(1)

# Initialize GMSH
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

# Create output directory
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

try:
    # Load the mesh file
    gmsh.open(mesh_file)
    logging.info(f"Successfully opened the mesh file: {mesh_file}")

    # Get physical groups
    physical_groups = gmsh.model.getPhysicalGroups()
    logging.info(f"Physical Groups: {physical_groups}")

    # Process each physical group
    for dim, tag in physical_groups:
        name = gmsh.model.getPhysicalName(dim, tag) or f"group_{dim}_{tag}"
        filename = os.path.join(output_dir, f"{name.lower()}.dat")

        # Get entities in the physical group
        entity_tags = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        logging.info(f"Physical Group {name} (dim={dim}, tag={tag}) has entities: {entity_tags}")

        with open(filename, "w") as file:
            file.write(f"# Physical Group: {name} (dim={dim}, tag={tag})\n")

            for entity_tag in entity_tags:
                try:
                    logging.info(f"Processing entity {entity_tag} in group {name}")

                    # Get nodes for the entity
                    node_tags, node_coords, _ = gmsh.model.mesh.getNodes(dim=dim, tag=entity_tag)
                    if len(node_tags) == 0:
                        logging.warning(f"No nodes found for entity {entity_tag}.")
                        continue

                    file.write(f"\n# Entity {entity_tag} Nodes:\n")
                    for i, tag in enumerate(node_tags):
                        coords = node_coords[3 * i: 3 * i + 3]
                        file.write(f"{tag}: {coords}\n")

                    # Get elements for the entity
                    element_types, element_tags, element_node_tags = gmsh.model.mesh.getElements(dim=dim, tag=entity_tag)
                    if len(element_types) == 0:
                        logging.warning(f"No elements found for entity {entity_tag}.")
                        continue

                    file.write(f"\n# Entity {entity_tag} Elements:\n")
                    for etype, etags, enodes in zip(element_types, element_tags, element_node_tags):
                        file.write(f"Type {etype}:\n")
                        num_nodes_per_element = gmsh.model.mesh.getElementProperties(etype)[3]
                        for e, nodes in zip(etags, zip(*[iter(enodes)] * num_nodes_per_element)):
                            file.write(f"{e}: {nodes}\n")
                except Exception as e:
                    logging.error(f"Error processing entity {entity_tag} in group {name}: {e}")

            logging.info(f"Saved data for {name} to {filename}")

except Exception as e:
    logging.error(f"An error occurred: {e}")

finally:
    # Finalize GMSH
    gmsh.finalize()
    logging.info("GMSH finalized.")
    exit(0)
