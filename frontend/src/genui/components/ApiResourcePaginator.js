import React from 'react';
import Pagination from 'react-js-pagination';
import withUnmounted from '@ishawnwang/withunmounted';

class ApiResourcePaginator extends React.Component {
  abort = new AbortController();
  hasUnmounted = false;

  constructor(props) {
    super(props);
    this.state = {
      activePage: 1,
      activePageItems: [],
      totalCount: 0,
    };
  }

  componentDidMount() {
    this.fetchPage(this.props.url, this.state.activePage);
  }

  handlePageChange(pageNumber) {
    this.fetchPage(this.props.url, pageNumber);
  }

  fetchPage = (rootUrl, pageNumber) => {
    const url = `${rootUrl.toString()}?page=${pageNumber}`;
    fetch(url, {signal : this.abort.signal})
      .then(response => response.json())
      .then(data => {
        if (this.hasUnmounted) {
          return
        }
        this.setState({
          activePageItems: data.results,
          totalCount: data.count,
          activePage: pageNumber
        })
      })
      .catch(e => console.log(e))
  };

  render() {
    return (
      <React.Fragment>
        {this.props.children(this.state.activePageItems)}
        <Pagination
          activePage={this.state.activePage}
          itemsCountPerPage={this.props.itemsPerPage ? this.props.itemsPerPage : 10}
          totalItemsCount={this.state.totalCount}
          pageRangeDisplayed={5}
          onChange={this.handlePageChange.bind(this)}
          itemClass="page-item"
          linkClass="page-link"
        />
      </React.Fragment>
    );
  }
}

export default withUnmounted(ApiResourcePaginator);