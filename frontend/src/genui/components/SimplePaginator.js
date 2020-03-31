import React from 'react';
import JwPagination from 'jw-react-pagination';
import { smoothScrollToTop } from '../utils';

class SimplePaginator extends React.Component {
  constructor(props) {
    super(props);

    this.onChangePage = this.onChangePage.bind(this);

    this.state = {
      items: [],
      pageOfItems: []
    };
  }

  componentDidMount() {
    const items = [...this.props.items];
    this.setState({
      items,
      pageOfItems: []
    })
  }

  componentDidUpdate(prevProps, prevState, snapshot) {
    if (this.props.forceUpdate || (this.props.items.length !== prevState.items.length)) {
      const items = [...this.props.items];
      this.setState({
        items,
        pageOfItems: []
      })
    }
  }

  onChangePage(pageOfItems) {
    this.setState({ pageOfItems: pageOfItems });
    smoothScrollToTop();
  }

  render() {
    const Component = this.props.component;
    return (
      <div id="scrollTop" className="simple-paginator">
        {this.state.pageOfItems.map(item =>
          <Component {...this.props} key={item.id} pageItem={item}/>
        )}
        <JwPagination items={this.state.items} onChangePage={this.onChangePage} />
      </div>
    );
  }
}

export default SimplePaginator;