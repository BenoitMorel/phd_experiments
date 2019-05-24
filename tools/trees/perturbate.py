from ete3 import Tree
import random
import sys

def get_random_node(tree):
    l = []
    for node in tree.traverse():
        l.append(node)
    rdm = random.randint(1,1000000)
    return l[rdm % len(l)]
    
def get_random_son(node):
    rdm = random.randint(1, 1000000)
    return node.get_children()[rdm % len(node.get_children())]

def perturbate_tree(tree):
    node = get_random_node(tree)
    if (node.is_root()):
        return
    parent_swap = node.get_sisters()[0]
    if (parent_swap.is_leaf()):
        return
    parent_node = node.get_ancestors()[0]
    swap = get_random_son(parent_swap) 
    node.detach()
    swap.detach()
    parent_swap.add_child(node)
    parent_node.add_child(swap)

def perturbate(input_tree_path, output_tree_path, min_aRF):
    original_tree = Tree(open(input_tree_path).read(), format=1)
    tree = Tree(open(input_tree_path).read(), format=1)
    for i in range(0, 10000):
        perturbate_tree(tree)
        RF = tree.robinson_foulds(original_tree, unrooted_trees=True)
        aRF = float(RF[0]) / float(RF[1])
        print("RF: " + str(RF[0]) + " " + str(aRF))
        if (aRF >= float(min_aRF)):
            break
    print("Saving perturbated tree in " + output_tree_path)
    open(output_tree_path, "w").write(tree.write())

if (__name__ == "__main__"):
    if (len(sys.argv) != 4):
        print("Syntax: input_tree_path output_tree_path min_aRF")
        exit(1)
    input_tree_path = sys.argv[1]
    output_tree_path = sys.argv[2]
    min_aRF = float(sys.argv[3])
    perturbate(input_tree_path, output_tree_path, min_aRF)
    
