"""Submission ZIP must contain all three tasks with consecutive ANGLE indices."""
from pathlib import Path
from zipfile import ZipFile
import pytest
from scripts.package_submission import create_submission

def test_zip_contents(tmp_path):
    samples = tmp_path / "sample_submission"
    df = samples / "dfthcsos"
    df.mkdir(parents=True)
    (samples / "dqi.kmb").write_bytes(b"ref DQI")
    (df / "submission.npz").write_bytes(b"ref factors")
    for i in range(3):
        (df / f"angle_k_{i}.kmb").write_bytes(bytes([i]))
    opt = tmp_path / "hybrid_8_26.kmb"
    opt.write_bytes(b"optimized invsqrt")
    z = tmp_path / "submission.zip"
    members = create_submission(opt, tmp_path, z)
    assert len(members) == 6
    with ZipFile(z) as archive:
        assert archive.read("invsqrt.kmb") == b"optimized invsqrt"
        assert archive.read("dqi.kmb") == b"ref DQI"
        assert all(f"dfthcsos/angle_k_{i}.kmb" in archive.namelist() for i in range(3))

def test_missing_reference_files_rejected(tmp_path):
    opt = tmp_path / "candidate.kmb"
    opt.write_bytes(b"x")
    with pytest.raises(FileNotFoundError):
        create_submission(opt, tmp_path, tmp_path / "bad.zip")

def test_reject_discontinuous_angles(tmp_path):
    samples=tmp_path/"sample_submission"
    df=samples/"dfthcsos"
    df.mkdir(parents=True)
    (samples/"dqi.kmb").write_bytes(b"a")
    (df/"submission.npz").write_bytes(b"a")
    (df/"angle_k_0.kmb").write_bytes(b"a")
    (df/"angle_k_2.kmb").write_bytes(b"a")
    opt=tmp_path/"candidate.kmb"
    opt.write_bytes(b"b")
    with pytest.raises(ValueError,match="Noncontiguous"):
        create_submission(opt,tmp_path,tmp_path/"bad.zip")
