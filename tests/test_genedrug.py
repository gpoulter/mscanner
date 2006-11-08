#!env python

from path import path
import unittest
from pprint import PrettyPrinter
from genedrug import *
from article import Article

class GeneDrugFilterTests(unittest.TestCase):
    """Tests for GeneDrugFilter class

    Tests: listSentences
    """
    def test_sentences(self):
        sentence_text="One A. Dr. Smith. Two B. A question? An exclamation! One. lowercase. Just... an ellipsis."
        sentences_correct=[
            ('One A. ', 0, 7),
            ('Dr. Smith. ', 7, 18),
            ('Two B. ', 18, 25),
            ('A question? ', 25, 37),
            ('An exclamation! ', 37, 53),
            ('One. lowercase. ', 53, 69),
            ('Just... an ellipsis.', 69, 89)]

if __name__ == "__main__":
    unittest.main()

