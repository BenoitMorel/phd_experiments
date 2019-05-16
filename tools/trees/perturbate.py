from ete3 import Tree
import random

def get_random_node(tree):
    l = []
    for node in tree.traverse():
        l.append(node)
    rdm = random.randint(1,1000000)
    return l[rdm % len(l)]
    
def perturbate(tree):
    
    node = get_random_node(tree)
    if (node.is_root()):
        return
    sister = node.get_sisters()[0]
    node.detach()
    sister.add_child(node)


species_tree = open("../BenoitDatasets/families/swiss/speciesTree.newick").read()
tree1 = Tree(species_tree, format=1) 
tree2 = Tree(species_tree, format=1) 

print(tree1.robinson_foulds(tree2, unrooted_trees = True)[0])
for i in range(0, 100):
    perturbate(tree1)
    print(tree1.robinson_foulds(tree2, unrooted_trees = True)[0])


