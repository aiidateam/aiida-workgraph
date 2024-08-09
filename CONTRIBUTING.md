
### Running tests

For running the tests you require only need start rabbitmq in the background.
The profile that is created during the tests will automatically check its configuration.
To run all tests use

```console
pytest
```

To run only backend tests run

```console
pytest -m backend
```

To run only frontend tests
```console
pytest -m frontend
```

To not run these tests you can use the markers in the following way

```console
pytest -m "not backend and not frontend"
```

To debug the frontend tests you often want to see what happens in the tests.
By default they are run in headless mode, so no browser is shown.
To run the frontend tests in headed mode for you have to set an environment variable like this
```console
PYTEST_PLAYWRIGHT_HEADLESS=no pytest -m frontend
```

### Development on the GUI

For the development on the GUI we use the [REACT](https://react.dev) library
which can automatically refresh on changes of the JS files. To start the backend
server please run

```console
python aiida_workgraph/web/backend/main.py
```

then start the frontend server with
```console
npm --prefix aiida_workgraph/web/frontend start
```

The frontend server will refresh

### Troubleshooting

#### Tests are not updating after changes in code

You might want to clean your cache

```console
npm --prefix aiida_workgraph/web/frontend cache clean
```

and also clear your browsers cache or try to start new private window.
