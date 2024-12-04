def test_load_config():
    from aiida_workgraph.config import load_config

    config = load_config()
    assert isinstance(config, dict)
    assert config == {}
