# coding: utf8

# Import modules
# import collections
# import cfdtools.meshbase._mesh as _mesh
# import cfdtools.meshbase._connectivity as _conn

# import os

# import cfdtools.api as api
# import numpy as np

# from operator import itemgetter

# Gmsh element types to canonical types see description below after ReaderGmsh class object
gmshtype_elem = {
    1: "bar2",
    2: "tri3",
    3: "quad4",
    4: "tetra4",
    5: "hexa8",
    6: "prism6",
    7: "pyra5",
    8: "bar3",
    9: "tri6",
    10: "quad9",
    11: "tetra10",
    12: "hexa27",
    13: "prism18",
    14: "pyra14",
    15: "node1",
}
# # Actual number of vertices for a given cell type
# nodes_per_cell = {
#     "bi": 2,
#     "tri": 3,
#     "qua": 4,
#     "tet": 4,
#     "hex": 8,
#     "pri": 6,
#     "pyr": 5,
# }

# etype_from_gmsh = {
#     "lin" : "bar2",
#     "tri" : "tri3",
#     "hex" : "hexa8"
# }

# #================================================================================================
# #================================================================================================

# From GMSH doc -
#                                                     (node associations)
# *  1 :   2-node 1-D            line.
# *  2 :   3-node 2-D            triangle.
# *  3 :   4-node 2-D            quadrangle.
# *  4 :   4-node 3-D            tetrahedron.
# *  5 :   8-node 3-D            hexahedron.
# *  6 :   6-node 3-D            prism.
# *  7 :   5-node 3-D            pyramid.
# *  8 :   3-node 1-D 2-nd order line                 (2 w/vertices,  1 w/edge).
# *  9 :   6-node 2-D 2-nd order triangle             (3 w/vertices,  3 w/edges).
# * 10 :   9-node 2-D 2-nd order quadrangle           (4 w/vertices,  4 w/edges,  1 w/face).
# * 11 :  10-node 3-D 2-nd order tetrahedron          (4 w/vertices,  6 w/edges).
# * 12 :  27-node 3-D 2-nd order hexahedron           (8 w/vertices, 12 w/edges,  6 w/faces,  1 w/volume).
# * 13 :  18-node 3-D 2-nd order prism                (6 w/vertices,  9 w/edges,  3 w/quad faces).
# * 14 :  14-node 3-D 2-nd order pyramid              (5 w/vertices,  8 w/edges,  1 w/quad face).
#   15 :   1-node 0-D            point.
# * 16 :   8-node 2-D 2-nd order quadrangle           (4 w/vertices,  4 w/edges).
# * 17 :  20-node 3-D 2-nd order hexahedron           (8 w/vertices, 12 w/edges).
# * 18 :  15-node 3-D 2-nd order prism                (6 w/vertices,  9 w/edges).
# * 19 :  13-node 3-D 2-nd order pyramid              (5 w/vertices,  8 w/edges).
# * 20 :   9-node 2-D 3-rd order incomplete triangle  (3 w/vertices,  6 w/edges)
# * 21 :  10-node 2-D 3-rd order triangle             (3 w/vertices,  6 w/edges,  1 w/face)
# * 22 :  12-node 2-D 4-th order incomplete triangle  (3 w/vertices,  9 w/edges)
# * 23 :  15-node 2-D 4-th order triangle             (3 w/vertices,  9 w/edges,  3 w/face)
#   24 :  15-node 2-D 5-th order incomplete triangle  (3 w/vertices, 12 w/edges)
#   25 :  21-node 2-D 5-th order complete triangle    (3 w/vertices, 12 w/edges,  6 w/face)
# * 26 :   4-node 1-D 3-rd order edge                 (2 w/vertices,  2 w/edge)
#   27 :   5-node 1-D 4-th order edge                 (2 w/vertices,  3 w/edge)
#   28 :   6-node 1-D 5-th order edge                 (2 w/vertices,  4 w/edge)
#   29 :  20-node 3-D 3-rd order tetrahedron          (4 w/vertices, 12 w/edges,  4 w/faces)
#   30 :  35-node 3-D 4-th order tetrahedron          (4 w/vertices, 18 w/edges, 12 w/faces,  1 w/volume)
#   31 :  56-node 3-D 5-th order tetrahedron          (4 w/vertices, 24 w/edges, 24 w/faces,  4 w/volume)
#   92 :  64-node 3-D 3-rd order hexahedron           (8 w/vertices, 24 w/edges, 24 w/faces,  8 w/volume)
#   93 : 125-node 3-D 4-th order hexahedron           (8 w/vertices, 36 w/edges, 54 w/faces, 27 w/volume)


# *  1:Line:              *  8:Line3:      * 26:Line4:

# 0----------1 --> u      0-----2----1     0----2----3----1


# *  2:Triangle:          *  9:Triangle6:      * 20/21:Triangle9/10:  * 22/23:Triangle12/15:

# v
# ^                                                                   2
# |                                                                   | \
# 2                       2                    2                      9   8
# |`\                     |`\                  | \                    |     \
# |  `\                   |  `\                7   6                 10 (14)  7
# |    `\                 5    `4              |     \                |         \
# |      `\               |      `\            8  (9)  5             11 (12) (13) 6
# |        `\             |        `\          |         \            |             \
# 0----------1 --> u      0-----3----1         0---3---4---1          0---3---4---5---1


# *  3:Quadrangle:       * 16:Quadrangle8:       * 10:Quadrangle9:

#       v
#       ^
#       |
# 3-----------2          3-----6-----2           3-----6-----2
# |     |     |          |           |           |           |
# |     |     |          |           |           |           |
# |     +---- | --> u    7           5           7     8     5
# |           |          |           |           |           |
# |           |          |           |           |           |
# 0-----------1          0-----4-----1           0-----4-----1


# *  4:Tetrahedron:                     * 11:Tetrahedron10:

#                    v
#                  .
#                ,/
#               /
#            2                                     2
#          ,/|`\                                 ,/|`\
#        ,/  |  `\                             ,/  |  `\
#      ,/    '.   `\                         ,6    '.   `5
#    ,/       |     `\                     ,/       8     `\
#  ,/         |       `\                 ,/         |       `\
# 0-----------'.--------1 --> u         0--------4--'.--------1
#  `\.         |      ,/                 `\.         |      ,/
#     `\.      |    ,/                      `\.      |    ,9
#        `\.   '. ,/                           `7.   '. ,/
#           `\. |/                                `\. |/
#              `3                                    `3
#                 `\.
#                    ` w


# *  5:Hexahedron:        * 17:Hexahedron20:     * 12:Hexahedron27:

#        v
# 3----------2            3----13----2           3----13----2
# |\     ^   |\           |\         |\          |\         |\
# | \    |   | \          | 15       | 14        |15    24  | 14
# |  \   |   |  \         9  \       11 \        9  \ 20    11 \
# |   7------+---6        |   7----19+---6       |   7----19+---6
# |   |  +-- |-- | -> u   |   |      |   |       |22 |  26  | 23|
# 0---+---\--1   |        0---+-8----1   |       0---+-8----1   |
#  \  |    \  \  |         \  17      \  18       \ 17    25 \  18
#   \ |     \  \ |         10 |        12|        10 |  21    12|
#    \|      w  \|           \|         \|          \|         \|
#     4----------5            4----16----5           4----16----5


# *  6:Prism:                 * 18:Prism15:          * 13:Prism18:

#            w
#            ^
#            |
#            3                       3                      3
#          ,/|`\                   ,/|`\                  ,/|`\
#        ,/  |  `\               12  |  13              12  |  13
#      ,/    |    `\           ,/    |    `\          ,/    |    `\
#     4------+------5         4------14-----5        4------14-----5
#     |      |      |         |      8      |        |      8      |
#     |    ,/|`\    |         |      |      |        |    ,/|`\    |
#     |  ,/  |  `\  |         |      |      |        |  15  |  16  |
#     |,/    |    `\|         |      |      |        |,/    |    `\|
#    ,|      |      |\        10     |      11       10-----17-----11
#  ,/ |      0      | `\      |      0      |        |      0      |
# u   |    ,/ `\    |    v    |    ,/ `\    |        |    ,/ `\    |
#     |  ,/     `\  |         |  ,6     `7  |        |  ,6     `7  |
#     |,/         `\|         |,/         `\|        |,/         `\|
#     1-------------2         1------9------2        1------9------2


# *  7:Pyramid:                * 19:Pyramid13:              * 14:Pyramid14:

#                4                            4                            4
#              ,/|\                         ,/|\                         ,/|\
#            ,/ .'|\                      ,/ .'|\                      ,/ .'|\
#          ,/   | | \                   ,/   | | \                   ,/   | | \
#        ,/    .' | `.                ,/    .' | `.                ,/    .' | `.
#      ,/      |  '.  \             ,7      |  12  \             ,7      |  12  \
#    ,/       .' w |   \          ,/       .'   |   \          ,/       .'   |   \
#  ,/         |  ^ |    \       ,/         9    |    11      ,/         9    |    11
# 0----------.'--|-3    `.     0--------6-.'----3    `.     0--------6-.'----3    `.
#  `\        |   |  `\    \      `\        |      `\    \     `\        |      `\    \
#    `\     .'   +----`\ - \ -> v  `5     .'        10   \      `5     .' 13     10   \
#      `\   |    `\     `\  \        `\   |           `\  \       `\   |           `\  \
#        `\.'      `\     `\`          `\.'             `\`         `\.'             `\`
#           1----------------2            1--------8-------2           1--------8-------2
#                     `\
#                        u

###################################################################################################
