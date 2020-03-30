import React from 'react';

class LogInManager extends React.Component {
  // FIXME: make sure we don't update state when this component is unmounted

  constructor(props) {
    super(props);

    this.state = {
      errors: [],
      loginSuccess: false
    }
  }

  componentDidMount() {
    if (this.props.fetchUser) {
      this.props.fetchUserInfo((userData) => {
        if (userData) {
          this.props.setUser(userData);
        }
      });
    }
  }

  sendLogInRequest = (credentials) => {
    fetch(new URL('login/', this.props.apiUrls.accountsRoot), {
      credentials: "include",
      "headers": {
        "Accept": "application/json",
        "Content-Type": "application/json",
      },
      "body": JSON.stringify(credentials),
      "method": "POST"
    }).then(response => response.json())
      .then(data => {
        if (data.hasOwnProperty('non_field_errors')) {
          this.setState(prevState => ({errors: prevState.errors.concat(data['non_field_errors'])}));
          return null;
        }
        return data;
      })
      .then(data => {
        if(data) {
          this.setState({loginSuccess: true});
          this.props.fetchUserInfo((userData) => {
            this.props.setUser(userData);
          });
        }
      })
      .catch(e => console.log(e));
  };

  sendLogOutRequest = () => {
    fetch(new URL('logout/', this.props.apiUrls.accountsRoot), {
      "credentials": "include",
      "headers": {
        "Accept": "application/json",
        "Content-Type": "application/json",
      },
      "method": "POST"
    }).then(response => response.json())
      .then(data => {
        // console.log(data);
        this.setState({loginSuccess: false});
        this.props.setUser(null)
      })
      .catch(e => console.log(e))
  };

  render() {
    if (this.props.component) {
      const LogInForm = this.props.component;
      return (
        <LogInForm
          {...this.props}
          onSubmit={this.sendLogInRequest}
          logoutUser={this.sendLogOutRequest}
          apiErrors={this.state.errors}
          loginSuccess={this.state.loginSuccess}
        />
      )
    } else {
      return this.props.children(this.sendLogInRequest, this.sendLogOutRequest, this.state.loginSuccess, this.state.errors)
    }
  }
}

export default LogInManager;