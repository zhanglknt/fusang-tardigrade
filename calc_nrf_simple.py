#!/usr/bin/env python3
"""Pure Python nRF (normalized Robinson-Foulds) calculator. No dependencies.

Parses Newick trees, extracts bipartitions, computes nRF distance.
"""

import sys
import re


def parse_newick(nwk_str):
    """Parse Newick tree, return set of leaf names and list of bipartitions.
    
    A bipartition is a frozenset of leaf names on one side of an internal edge.
    """
    # Remove comments in square brackets
    nwk_str = re.sub(r'\[[^\]]*\]', '', nwk_str)
    # Remove trailing semicolon
    nwk_str = nwk_str.strip().rstrip(';')
    
    # Tokenize: split on ( ) , but keep delimiters
    # Also handle : for branch lengths
    tokens = []
    i = 0
    while i < len(nwk_str):
        c = nwk_str[i]
        if c in '(),':
            tokens.append(c)
            i += 1
        elif c in ' \t\n\r':
            i += 1  # skip whitespace
        elif c == ':':
            # Branch length - read until , ) or end
            i += 1
            start = i
            while i < len(nwk_str) and nwk_str[i] not in '(),':
                i += 1
            # Skip the branch length value, just note there was one
            # tokens.append((':', nwk_str[start:i]))
        else:
            # Leaf name or internal node label
            start = i
            while i < len(nwk_str) and nwk_str[i] not in '(),:\n\r\t ':
                i += 1
            tokens.append(nwk_str[start:i])
    
    # Parse with stack
    # Each stack entry: (children_names, node_label)
    # children_names: list of sets of leaf names reachable from each child
    stack = []
    bipartitions = set()
    all_leaves = set()
    
    for token in tokens:
        if token == '(':
            stack.append([])  # Start a new subtree, collect children
        elif token == ',':
            pass  # Separator between siblings
        elif token == ')':
            # This subtree is complete
            # The attached label (if any) will be the next token
            pass
        else:
            # A label: could be leaf or internal node
            label = token
            
            if not stack:
                # Single node tree (edge case)
                all_leaves.add(label)
                continue
            
            if stack and isinstance(stack[-1], list):
                # This is a child of the current open group
                current_children = stack[-1]
                
                # Check if this is a leaf or internal node
                # If the previous token was ')', this label belongs to that internal node
                # Otherwise it's a leaf
                
                # Actually, in our tokenization, a leaf is just a label
                # For now, treat as leaf
                leaf_set = {label}
                all_leaves.add(label)
                current_children.append(leaf_set)
    
    # Now post-process: collapse the tree and extract bipartitions
    def resolve_subtree(children_list):
        """Resolve a list of child leaf sets into a set of leaf sets for each subtree.
        Returns: list of (child_leaf_set, bipartition_set_to_add) for each child.
        """
        if not children_list:
            return set()
        
        # Each child is currently a set of leaves
        # For each internal edge between a child and the parent,
        # the bipartition is (child_leaf_set, all_other_leaves - child_leaf_set)
        result_leaves = set()
        for child_leaves in children_list:
            result_leaves |= child_leaves
        
        return result_leaves
    
    # This approach is getting too complicated. Let me use a recursive descent parser instead.
    return all_leaves, bipartitions


def _parse_subtree(s, pos):
    """Recursive descent Newick parser.
    Returns: (leaf_set_for_subtree, bipartitions_list, new_pos)
    """
    bipartitions = []
    leaves = set()
    
    if pos >= len(s):
        return leaves, bipartitions, pos
    
    if s[pos] == '(':
        # Internal node
        children_leaves = []
        pos += 1
        
        while pos < len(s):
            child_leaves, child_bps, pos = _parse_subtree(s, pos)
            children_leaves.append(child_leaves)
            bipartitions.extend(child_bps)
            
            if pos >= len(s):
                break
            if s[pos] == ')':
                pos += 1
                break
            elif s[pos] == ',':
                pos += 1
            else:
                break
        
        # Combine all children leaves
        leaves = set()
        for cl in children_leaves:
            leaves |= cl
        
        # Skip label (if any)
        while pos < len(s) and s[pos] not in '(),':
            pos += 1
        # Skip branch length (if any)
        if pos < len(s) and s[pos] == ':':
            pos += 1
            while pos < len(s) and s[pos] not in '(),':
                pos += 1
        
        # Generate bipartitions for each edge from this internal node to children
        all_leaves_set = leaves
        for child_leaves in children_leaves:
            if child_leaves and child_leaves != all_leaves_set:
                side = frozenset(child_leaves)
                other = frozenset(all_leaves_set - child_leaves)
                # Canonical: always store the side with fewer leaves or lexicographically smaller
                if len(side) < len(other) or (len(side) == len(other) and sorted(side) < sorted(other)):
                    bipartitions.append(side)
                else:
                    bipartitions.append(other)
        
        return leaves, bipartitions, pos
    
    else:
        # Leaf node
        name_start = pos
        while pos < len(s) and s[pos] not in '(),:':
            pos += 1
        name = s[name_start:pos].strip()
        
        # Skip branch length
        if pos < len(s) and s[pos] == ':':
            pos += 1
            while pos < len(s) and s[pos] not in '(),':
                pos += 1
        
        return {name} if name else set(), [], pos


def get_bipartitions_from_newick(nwk_str):
    """Extract all bipartitions from a Newick string."""
    # Clean and normalize
    nwk_str = re.sub(r'\[[^\]]*\]', '', nwk_str)
    nwk_str = nwk_str.strip().rstrip(';')
    
    leaves, bipartitions, pos = _parse_subtree(nwk_str, 0)
    
    # Remove empty bipartitions
    bipartitions = [bp for bp in bipartitions if len(bp) > 0]
    
    return leaves, set(bipartitions)


def calc_nrf(tree1_path, tree2_path):
    """Calculate normalized Robinson-Foulds distance."""
    with open(tree1_path, encoding='utf-8', errors='ignore') as f:
        nwk1 = f.read()
    with open(tree2_path, encoding='utf-8', errors='ignore') as f:
        nwk2 = f.read()
    
    leaves1, bp1 = get_bipartitions_from_newick(nwk1)
    leaves2, bp2 = get_bipartitions_from_newick(nwk2)
    
    all_leaves = leaves1 | leaves2
    
    # Verify same leaf set
    if leaves1 != leaves2:
        try:
            print(f"WARNING: Leaf sets differ! |tree1|={len(leaves1)}, |tree2|={len(leaves2)}")
            only1_list = [str(x) for x in sorted(list(leaves1 - leaves2))[:5]]
            only2_list = [str(x) for x in sorted(list(leaves2 - leaves1))[:5]]
            print(f"  Only in tree1: {only1_list}{'...' if len(leaves1 - leaves2) > 5 else ''}")
            print(f"  Only in tree2: {only2_list}{'...' if len(leaves2 - leaves1) > 5 else ''}")
        except UnicodeEncodeError:
            print(f"WARNING: Leaf sets differ! |tree1|={len(leaves1)}, |tree2|={len(leaves2)}")
    
    shared = bp1 & bp2
    only1 = bp1 - bp2
    only2 = bp2 - bp1
    
    rf = len(only1) + len(only2)
    n_taxa = len(all_leaves)
    max_rf = 2 * (n_taxa - 3) if n_taxa > 3 else 0
    
    nrf = rf / max_rf if max_rf > 0 else 0.0
    
    print(f"Tree1: {tree1_path}")
    print(f"Tree2: {tree2_path}")
    print(f"Taxa: {n_taxa}")
    print(f"Bipartitions: tree1={len(bp1)}, tree2={len(bp2)}")
    print(f"Shared: {len(shared)}")
    print(f"RF distance: {rf} (= {len(only1)} + {len(only2)})")
    print(f"Max RF: {max_rf}")
    print(f"nRF distance: {nrf:.6f}")
    return nrf


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} tree1.nwk tree2.nwk")
        sys.exit(1)
    calc_nrf(sys.argv[1], sys.argv[2])
