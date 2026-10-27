/** @type {import('tailwindcss').Config} */
module.exports = {
    content: [
        "./src/**/*.jl"
    ],
    darkMode: 'class',
    theme: {
        extend: {
            fontFamily: {
                sans: ['-apple-system', 'BlinkMacSystemFont', 'SF Pro Display', 'SF Pro Text', 'Helvetica Neue', 'Inter', 'Segoe UI', 'Roboto', 'Arial', 'sans-serif'],
            }
        }
    },
    plugins: [],
}
