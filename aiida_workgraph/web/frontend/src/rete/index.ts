import { createEditor as createDefaultEditor } from './default'
import { createEditor as createCustomEditor } from './customization'

const factory = {
  'default': createDefaultEditor,
  'customization': createCustomEditor,
}
// eslint-disable-next-line no-restricted-globals, no-undef
const query = typeof location !== 'undefined' && new URLSearchParams(location.search)
const name = ((query && query.get('template')) || 'default') as keyof typeof factory

const createEditor = factory[name]

if (!createEditor) {
  throw new Error(`template with name ${name} not found`)
}

export {
  createEditor
}
