import React from 'react';
import withUnmounted from '@ishawnwang/withunmounted';

class ComponentWithPagedResources extends React.Component {
  abort = new AbortController();
  hasUnmounted = false;

  constructor(props) {
    super(props);

    this.state = this.initState(0);
  }

  initState = (revision) => {
    const definition = this.props.definition;
    const data = {};
    Object.keys(definition).forEach(key => {
      data[key] = {
        items : [],
        lastPage : null,
        nextPage : definition[key],
        finished : false,
      }
    });
    return {
      data : data,
      isUpdating: true,
      revision: revision,
    }
  };

  componentDidMount() {
    this.updateAll();
  }

  updateAll = () => {
    this.setState({
      revision: this.state.revision + 1,
      isUpdating: true
    });

    const data = this.state.data;
    Object.keys(data).forEach(key => {
      this.fetchData(key, data[key].nextPage);
    });
  };

  // componentDidUpdate(prevProps, prevState, snapshot) {
  //   if (this.props.updateCondition) {
  //     if (this.props.updateCondition(prevProps, this.props, prevState, this.state, snapshot)) {
  //       if (this.hasUnmounted) {
  //         return
  //       }
  //       this.setState(
  //         prevState => this.initState(prevState.revision)
  //         , () => {
  //           this.updateAll();
  //         }
  //       );
  //     }
  //   }
  // }

  componentWillUnmount() {
    this.abort.abort();
  }

  isFinished(state) {
    let finished = true;
    Object.keys(state.data).forEach(key => {
      finished = finished && state.data[key].finished;
    });

    return finished;
  }

  fetchData = (key, page) => {
    fetch(page, {signal : this.abort.signal, credentials: "include"})
      .then(response => response.json())
      .then(data => {
        let nextPage = null;
        let finished = false;
        if (data.next) {
          nextPage = new URL(data.next);
        } else {
          finished = true;
        }

        if (this.hasUnmounted) {
          return
        }
        // console.log(page.toString());
        this.setState(prevState => {
          prevState.data[key].items = prevState.data[key].items.concat(data.results);
          prevState.data[key].lastPage = page;
          prevState.data[key].nextPage = nextPage;
          prevState.data[key].finished = finished;

          if (this.isFinished(prevState)) {
            prevState.isUpdating = false;
          }

          return prevState;
        }, () => {
          if (nextPage) {
            this.fetchData(key, nextPage);
          }
        })
      })
      .catch(e => {
        // console.log(e) FIXME: catch some important errors and only let abort signals and stuff pass
      })
  };

  render() {
    const ret = {};
    const data = this.state.data;
    Object.keys(data).forEach(key => {
      ret[key] = data[key].items;
    });
    return this.props.children(ret, !this.state.isUpdating, this.state.revision);
  }
}

export default withUnmounted(ComponentWithPagedResources);