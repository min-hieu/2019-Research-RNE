document.addEventListener('DOMContentLoaded', () => {

    const form = document.querySelector('form')
    document.querySelector('form').addEventListener('submit', event => {
      event.preventDefault()
      console.log('submit')
    })
    
    const inputField = document.querySelector('input[type="text"]')
    
    inputField.addEventListener('input', event => {
      console.log(`The value entered is ${inputField.value}`)
    })
    
    inputField.addEventListener('cut', event => {
      console.log('cut')
    })
    
    inputField.addEventListener('copy', event => {
      console.log('copy')
    })
    
    inputField.addEventListener('paste', event => {
      console.log('paste')
    })
  })