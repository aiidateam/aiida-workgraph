import pytest

import os


@pytest.fixture(scope="module")
def set_backend_server_settings(aiida_profile):
    os.environ["AIIDA_WORKGRAPH_GUI_PROFILE"] = aiida_profile.name
