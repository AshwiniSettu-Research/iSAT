from flask import Flask, render_template

# initialize app
app = Flask(__name__)

@app.route('/')
def main():
    return render_template('main.html')

@app.route('/index')
def index():
    return render_template('index.html')


if __name__ == '__main__':
    app.run(debug=True)
