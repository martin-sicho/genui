import React from "react";
import withUnmounted from '@ishawnwang/withunmounted';

class LiveObject extends React.Component {
  abort = new AbortController();
  hasUnmounted = false;

  constructor(props) {
    super(props);

    this.objURL = this.props.url;
    this.interval = this.props.updateInterval ? this.props.updateInterval: 2000;
    this.intervalID = null;
    this.state = {
      instance: null,
      isUpdating: true
    }
  }

  componentDidMount() {
    this.fetchUpdates();
  }

  componentWillUnmount() {
    this.abort.abort();
    clearTimeout(this.intervalID);
  }

  fetchUpdates = () => {
    fetch(this.objURL, {signal : this.abort.signal, credentials: "include",})
      .then(response => this.props.handleResponseErrors(response))
      .then((data) => {
        if (this.hasUnmounted) {
          return
        }
        this.setState({
          instance : data,
          isUpdating: false
        });
        this.intervalID = setTimeout(this.fetchUpdates, this.interval);
      })
      .catch(
        (error) => {
          console.log(error);
        }
      )
  };

  updateInstance = (data) => {
    if (this.hasUnmounted) {
      return
    }

    this.setState({
      isUpdating : true
    });
    const error_msg = `Failed to update resource: ${this.objURL}`;
    fetch(
      this.objURL
      , {
        method: 'PATCH'
        , signal : this.abort.signal
        , body: JSON.stringify(data)
        , headers: {
          'Content-Type': 'application/json'
        },
        credentials: "include",
      }
    )
      .then(response => this.props.handleResponseErrors(response, error_msg))
      .then(
        data => {
          if (this.hasUnmounted) {
            return
          }
          this.setState({
            instance : data,
            isUpdating : false
          })
        }
      )
      .catch(
        (error) => console.log(error)
      );
  };

  render() {

    if (this.state.instance === null) {
      return <div>Loading...</div>
    }

    return this.props.children(this.state.instance, this.updateInstance, this.state.isUpdating)
  }
}

// LiveObject.propTypes = {
//   children: PropTypes.func.isRequired
// };

export default withUnmounted(LiveObject);