#!/usr/bin/env python

"""
Created by: Lee Bergstrand (2024)

Description: Code for dealing with genome annotations.
"""


class AnnotationMap(object):
    """
    A class representing a map of the genome annotation bioinformatics software to be run on the genome.
    """

    def __init__(self, annotations):
        """
        Initialize the object with the provided annotations.
        Args:
            annotations (list[str]): A list annotation software to be run.
        Dynamic Attributes:
            dfast_func (bool): True if 'dfast_func' annotation should be run, False otherwise.
            eggnog (bool): True if 'eggnog' annotation should be run, False otherwise.
            gtdbtk (bool): True if 'gtdbtk' annotation should be run, False otherwise.
            checkm2 (bool): True if 'checkm2' annotation should be run, False otherwise.
            coverage (bool): True if 'coverage' annotation should be run, False otherwise.
        """
        self.expected_annotations = ['dfast_func', 'eggnog', 'gtdbtk', 'checkm2', 'coverage']

        # Convert the annotations to be lower cases in case the user makes a small typo.
        annotations_lower = [anno.lower() for anno in annotations]

        for current_annotation in self.expected_annotations:
            if current_annotation in annotations_lower:
                setattr(self, current_annotation, True)  # Set the attributes dynamically.
            else:
                setattr(self, current_annotation, False)

    def __repr__(self):
        """
        Generate a representation of the object state in the form of AnnotationMap(dfast_func=..., eggnog=..., ...).
        """
        repr_parts = []
        for current_annotation in self.expected_annotations:
            repr_parts.append(f"{current_annotation}={getattr(self, current_annotation)}")
        return f"AnnotationMap({', '.join(repr_parts)})"

    def __str__(self):
        return self.__repr__()
