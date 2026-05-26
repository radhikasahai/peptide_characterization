"""Tests for the FastAPI web backend."""

from __future__ import annotations

import pytest

fastapi = pytest.importorskip("fastapi")
httpx = pytest.importorskip("httpx")

from fastapi.testclient import TestClient

from api.main import app

client = TestClient(app)


def test_health():
    response = client.get("/health")
    assert response.status_code == 200
    assert response.json()["status"] == "ok"


def test_benchmark():
    response = client.get("/api/benchmark")
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["valid"] >= 1
    assert len(data["examples"]) >= 1


def test_characterize_gly_dipeptide():
    response = client.post(
        "/api/characterize",
        json={"sequence": "GG", "hydrolyze_sidechains": True},
    )
    assert response.status_code == 200
    data = response.json()
    assert data["valid"] is True
    assert data["smiles"] == "NCC(=O)NCC(=O)O"
    assert data["structure_png_base64"]
    assert data["descriptors"]["molecular_weight"] > 0


def test_characterize_invalid():
    response = client.post(
        "/api/characterize",
        json={"sequence": "AX", "hydrolyze_sidechains": True},
    )
    assert response.status_code == 400


def test_generate_motif():
    response = client.post(
        "/api/generate",
        json={"mode": "motif", "motif": "RGD", "repeats": 2, "n_term_flank": "A"},
    )
    assert response.status_code == 200
    assert response.json()["sequence"] == "ARGDRGD"
