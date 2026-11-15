"""Deterministic primitives shared by future distributed connectivity code."""

from dataclasses import dataclass
import hashlib
import struct


def canonical_face_key(face_type, global_node_ids):
    """Return an orientation-independent key for a face with global node IDs."""

    nodes = tuple(sorted(int(node) for node in global_node_ids))
    if not nodes:
        raise ValueError("a face must contain at least one node")
    if len(set(nodes)) != len(nodes):
        raise ValueError("a face key cannot contain duplicate node IDs")
    return str(face_type), nodes


def face_owner(face_key, comm_size):
    """Choose a deterministic owner rank without relying on Python hash()."""

    if comm_size < 1:
        raise ValueError("comm_size must be positive")
    face_type, nodes = face_key
    digest = hashlib.blake2b(digest_size=8)
    digest.update(face_type.encode("utf-8"))
    for node in nodes:
        digest.update(struct.pack("<q", node))
    return int.from_bytes(digest.digest(), byteorder="little") % comm_size


@dataclass(frozen=True)
class FaceCandidate:
    """One local cell face sent to the rank that owns its canonical key.

    The local orientation is preserved in local_node_ids; the key itself is
    orientation-independent and is used only to find matching candidates.
    """

    face_type: str
    global_node_ids: tuple
    global_cell_id: int
    local_cell_index: int
    local_node_ids: tuple

    @property
    def key(self):
        return canonical_face_key(self.face_type, self.global_node_ids)
