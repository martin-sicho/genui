import React from 'react';
import withUnmounted from '@ishawnwang/withunmounted';

class ComponentWithPagedResources extends React.Component {
  abort = new AbortController();
  hasUnmounted = false;

  constructor(props) {
    super(props);

    this.state = this.initState({});
  }

  initState = (prevState) => {
    const definition = this.props.definition;
    const data = {};
    Object.keys(definition).forEach(key => {
      data[key] = {
        items : [],
        lastPage : null,
        nextPage : definition[key],
      }
    });
    return {
      data : data
    }
  };

  componentDidMount() {
    this.updateAll();
  }

  updateAll = () => {
    const data = this.state.data;
    Object.keys(data).forEach(key => {
      this.fetchData(key, data[key].nextPage);
    });
  };

  componentDidUpdate(prevProps, prevState, snapshot) {
    if (this.props.updateCondition) {
      if (this.props.updateCondition(prevProps, this.props, this.state, prevState, snapshot)) {
        this.setState(
          prevState => this.initState(prevState)
          , () => {
            this.updateAll();
          }
        );
      }
    }
  }

  componentWillUnmount() {
    this.abort.abort();
  }

  fetchData = (key, page) => {
    fetch(page, {signal : this.abort.signal})
      .then(response => response.json())
      .then(data => {
        let nextPage = null;
        if (data.next) {
          nextPage = new URL(data.next);
        }

        if (this.hasUnmounted) {
          return
        }
        this.setState(prevState => {
          let items = prevState.data[key].items;
          items = items.concat(data.results);

          prevState.data[key] = {
            items : items,
            lastPage: page,
            nextPage: nextPage,
          };
          return prevState;
        }, () => {
          if (nextPage) {
            this.fetchData(key, nextPage);
          }
        })
      })
      .catch(e => console.log(e))
  };

  render() {
    const ret = {};
    const data = this.state.data;
    Object.keys(data).forEach(key => {
      ret[key] = data[key].items;
    });
    // console.log(ret);
    return this.props.children(ret);
  }
}

export default withUnmounted(ComponentWithPagedResources);