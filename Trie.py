# Inspired from https://albertauyeung.github.io/2020/06/15/python-trie.html/
# Hamming distance from https://github.com/volpato30/hamming-d-search/blob/master/trie.py


class TrieNode:
    """Node objects in Trie"""

    def __init__(self, char: str) -> None:
        self.char = char
        self.is_end = False
        self.counter = 0  # Counts number of times inserted, maybe want to change to something else?
        self.children = {}
        # Want to store anything else?
        # Want to insert some check on char?


class Trie:
    """Trie object storing TrieNodes"""

    def __init__(self) -> None:
        self.root = TrieNode("")

    def insert(self, word: str) -> None:
        node = self.root

        # Check if there is a child containing the character
        # If not, create new child with character
        for char in word:
            if char in node.children:
                node = node.children[char]
            else:
                # Create new node
                new_node = TrieNode(char)
                node.children[char] = new_node
                node = new_node  # continue with rest of the word

        # Entire word is inserted, mark that this is end of word, and increase counter
        node.is_end = True
        node.counter += 1

    def _dfs(self, node: TrieNode, prefix: str):
        """Depth-first traversal of Trie

        Args:
            - node: Node to start at
            - prefix: current prefix, stores word while traversing
        """

        # The end of a word has been reached
        # Add word and counter to output list
        if node.is_end:
            self.output.append((prefix + node.char, node.counter))

        for child in node.children.values():
            self.dfs(child, prefix + node.char)

    def query(self, x: str) -> list[tuple[str, int]]:
        """Finds all words that starts with prefix (x)"""

        # Create new empty list each time function is called
        self.output = []

        node = self.root

        # Checks that the word actually exists in the trie
        for char in x:
            if char in node.children:
                node = node.children[char]
            else:
                return []
        # Do depth first search, starting at last node in prefix (x)
        self.dfs(node, x[:-1])

        # return all words starting with prefix (x), ordered by least number of occurrence
        return sorted(self.output, key=lambda x: x[1])

    def search_hamming_dist(self, query_seq: str, d:int):
        result = []
        current_index = 0
        current_cost = 0

        node = self.root

        for child in node.children:
            self._search_recursive(node.children[child], child, query_seq, current_index, current_cost, d, result)

        print(result, current_index, current_cost)

    def _search_recursive(self, current_node: TrieNode, base: str, query_seq: str,
                          current_index: int, current_cost: int, d, result):
        if not (base == query_seq[current_index]):
            current_cost += 1
        if current_cost > d:
            return
        current_index += 1

        # stopping criteria
        if current_node.is_end:
            result.append(1) # TODO: Fix to add proper word or something interesting. Maybe add prefix to parameters
            return

        for base in current_node.children:
            self._search_recursive(current_node.children[base], base, query_seq, current_index, current_cost, d, result)



t = Trie()
t.insert("hey")
t.insert("hiy")
t.insert("yoo")

print(t.root.char)
t.search_hamming_dist("joy", 5)
